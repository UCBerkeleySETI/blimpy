r"""
Make waterfall plots of a large file.
"""

import os
import sys
from os.path import dirname, abspath, isdir
import gc
from argparse import ArgumentParser
from PIL import Image
import logging
logger_name = "stix"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.INFO)

# Plotting packages import
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("agg")

# Blimpy imports
from blimpy import Waterfall
from blimpy.plotting import plot_waterfall


def image_stitch(orientation, chunk_count, png_collection, path_saved_png):
    r""" Stitch together multiple PNGs into one
    
    Parameters
    ----------
    orientation : str
        Assembling images horizontally (h) or vertically (v)?
    chunk_count : int
        Number of chunks in the file.
    png_collection : list
        The set of PNG file paths whose images are to be stitched together.
    path_saved_png : str
        The path of where to save the final PNG file.
    """

    logger.info("Stitching together the images, orientation is {}".format(orientation))

    # Empty array to store loaded images.
    img_collection = []

    # Loop to load images.
    for ii in range(chunk_count):
        img = Image.open(png_collection[ii])
        if ii == 0:
            base_width, base_height = img.size
            logger.info("Base (first) image file pixel width = {}, height = {}"
                        .format(base_width, base_height))
        width, height = img.size
        if width != base_width or height != base_height:
            logger.warning("Skipping {} because pixel width = {}, height = {}, doesn't match the base!"
                           .format(png_collection[ii], width, height))
            continue
        img_collection.append(img)

    # Get dimensions of the last image (same as the rest).

    # Creating the template for the stitched image.
    if orientation == "h":
        new_img = Image.new("RGBA", (chunk_count * width, height))
    else:
        new_img = Image.new("RGBA", (width, chunk_count * height))

    # Initalising the template before pasting the images.
    init_width = 0
    init_height = 0

    # Pasting the images in sequence.
    counter = 0
    for img in img_collection:
        box = (init_width, init_height, init_width + width, init_height + height)
        try:
            new_img.paste(img, box)
        except ValueError:
            logger.error("stix: *** Oops, Image paste error, init_width={}, init_height={}, width={}, height={}"
                         .format(init_width, init_height, width, height))
            sys.exit(86)
        counter += 1
        if orientation == "h":
            init_width += width
        else:
            init_height += height

    # Save the output, re-sized.
    if orientation == "h":
        new_img = new_img.resize((chunk_count * width, height))
    else:
        new_img = new_img.resize((width, chunk_count * height))
    new_img.save(path_saved_png)
    logger.info("Concatenated {} images to {}".format(counter, path_saved_png))


def sort2(x, y):
    r""" Return lowest value, highest value"""
    if y < x:
        return y, x
    return x, y


def make_waterfall_plots(input_file, chunk_count, plot_dir, width, height, dpi, source_name=None):
    r"""
    Make waterfall plots of a given huge-ish file.
​
    Parameters
    ----------
    input_file : str
        Path of Filterbank or HDF5 input file to plot in a stacked mode.
    chunk_count : int
        The number of chunks to divide the entire bandwidth into.
    plot_dir : str
        Directory for storing the PNG files.
    width : float
        Plot width in inches.
    height : float
        Plot height in inches.
    dpi : int
        Plot dots per inch.
    source_name : str
        Source name from the file header.
    """

    # Get directory path for storing PNG file
    if plot_dir is None:
        dirpath = dirname(abspath(input_file)) + "/"
    else:
        if not isdir(plot_dir):
            os.mkdir(plot_dir)
        dirpath = plot_dir

    # Calculate frequency boundary variables.
    wf1 = Waterfall(input_file, load_data=False)
    nf = wf1.n_channels_in_file
    if nf % chunk_count != 0:
        msg = "Number of frequency chunks ({}) does not evenly divide the number of frequencies ({})!" \
              .format(chunk_count, nf)
        logger.warning(msg)
    fchunk_f_start = wf1.file_header["fch1"]
    fchunk_size = int(nf / chunk_count)
    fchunk_f_offset = fchunk_size * wf1.file_header["foff"]
    fchunk_f_stop = fchunk_f_start + (fchunk_size - 1) * wf1.file_header["foff"]

    # Produce the source_name to be used in the image title.
    if source_name is None:
        source_name = wf1.header["source_name"]
        if source_name.upper() == "UNKNOWN":
            source_name = wf1.header["rawdatafile"].replace(".0000.raw", "")

    # Begin PNG file collection.
    png_collection = []

    # Generate a plot/PNG for each frequency chunk.
    logger.info("Image width = {}, height = {}, dpi = {}"
                .format(width, height, dpi))
    for ii in range(0, chunk_count):

        ii_lowest, ii_highest = sort2(fchunk_f_start, fchunk_f_stop)

        # read in data
        wf = Waterfall(input_file, f_start=ii_lowest, f_stop=ii_highest)

        # Validate frequency range.
        logger.debug("stix: Processing chunk {} of {}, frequency lowest={}, highest={}"
                    .format(ii, chunk_count - 1, ii_lowest, ii_highest))

        # Make plot for this frequency chunk with plot_waterfall().
        wf.header["source_name"] = source_name + " chunk " + str(ii + 1) + " of " + str(chunk_count)
        plt.figure(wf.header["source_name"], figsize=(width, height), dpi=dpi)
        plot_waterfall(wf, f_start=ii_lowest, f_stop=ii_highest)

        # Save the figures.
        path_png = dirpath + source_name + "_chunk_{}".format(ii + 1) + ".png"
        png_collection.append(path_png)
        plt.savefig(path_png, dpi=dpi)
        logger.info("Saved plot: {}".format(path_png))

        # Delete Waterfall object.
        del wf
        gc.collect()
        plt.close("all")

        # Prepare for next chunk.
        fchunk_f_start += fchunk_f_offset
        fchunk_f_stop += fchunk_f_offset

    return png_collection, source_name


def cmd_tool(args=None):
    r"""Coomand line parser"""
    parser = ArgumentParser(description="Make waterfall plots from a single file.")
    parser.add_argument("input_file", type=str, default=None, help="Input file to plot")
    parser.add_argument("chunk_count", type=int, default=0, help="Count of same-sized frequency chunks")
    parser.add_argument("--plot_dir", "-p", type=str, default=".", help="Directory to receive the plot (.png).  Default: current directory.")
    parser.add_argument("--source_name", "-n", type=str, default=None, help="Source name to override the header field.  Default: None.")
    parser.add_argument("--stitch", "-s", type=str, default="n", choices=["n", "v", "h"],
                        help="Stitch files: n(no stitching), v(vertical), or h(horizontal).  Default: n.")
    parser.add_argument("--dpi", "-d", type=int, default=200,
                        help="Single image file dots per inch.  Default: 200.")
    parser.add_argument("--width", "-w", type=float, default=10,
                        help="Single image file image width in inches.  Default: 10.0.")
    parser.add_argument("--height", "-t", type=float, default=8,
                        help="Single image file image height in inches.  Default: 8.0.")

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    if args.input_file is None or args.chunk_count < 1:
        os.system("stix -h")
        sys.exit(0)
    if args.dpi < 50:
        logger.error("stix: --dpi must be at least 50 but I saw: {}!".format(args.dpi))
        os.system("stix -h")
        sys.exit(86)
    if args.width < 6:
        logger.error("stix: --width must be at least 6 but I saw: {}!".format(args.width))
        os.system("stix -h")
        sys.exit(86)
    if args.height < 5:
        logger.error("stix: --height must be at least 5 but I saw: {}!".format(args.height))
        os.system("stix -h")
        sys.exit(86)

    # Make the plots.
    logger.info("Begin")
    args.plot_dir += "/"
    args.input_file = os.path.abspath(args.input_file)
    png_collection, source_name = make_waterfall_plots(args.input_file,
                                                       args.chunk_count,
                                                       args.plot_dir,
                                                       args.width,
                                                       args.height,
                                                       args.dpi,
                                                       source_name=args.source_name)
    if args.stitch != "n":
        path_saved_png = args.plot_dir + source_name + ".png"
        image_stitch(args.stitch, args.chunk_count, png_collection, path_saved_png)
        for file in png_collection:
            os.remove(file)

    logger.info("End")


if __name__ == "__main__":
    cmd_tool()
