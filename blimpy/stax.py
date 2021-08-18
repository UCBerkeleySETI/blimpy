r"""
Make waterfall plots of s set of files.
"""

import os
import sys
from os.path import isdir
from argparse import ArgumentParser
import logging
logger_name = "stax"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.INFO)
from PIL import Image

# Plotting packages import
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("agg")

# Blimpy imports
from blimpy import Waterfall
from blimpy.plotting import plot_waterfall


def image_stitch(orientation, png_collection, path_saved_png):
    r""" Stitch together multiple PNGs into one"""

    logger.info("Stitching together the images, orientation is {}".format(orientation))

    # Empty array to store loaded images.
    img_collection = []
    num_pngs = len(png_collection)

    # Loop to load images.
    ii = 0
    for png_file in png_collection:
        img = Image.open(png_file)
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
        ii += 1

    # Get dimensions of the last image (same as the rest).

    # Creating the template for the stitched image.
    if orientation == "h":
        new_img = Image.new("RGBA", (num_pngs * width, height))
    else:
        new_img = Image.new("RGBA", (width, num_pngs * height))

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
            logger.error("stax: *** Oops, Image paste error, init_width={}, init_height={}, width={}, height={}"
                         .format(init_width, init_height, width, height))
            sys.exit(86)
        counter += 1
        if orientation == "h":
            init_width += width
        else:
            init_height += height

    # Save the output, re-sized.
    if orientation == "h":
        new_img = new_img.resize((num_pngs * width, height))
    else:
        new_img = new_img.resize((width, num_pngs * height))
    new_img.save(path_saved_png)
    logger.info("Concatenated {} images to {}".format(counter, path_saved_png))


def make_waterfall_plots(file_list, plot_dir, width, height, dpi, f_start, f_stop):
    r"""
    Make waterfall plots of a given huge-ish file.
​
    Parameters
    ----------
    file_list : list
        List of paths of Filterbank or HDF5 input file to plot.
    plot_dir : str
        Directory for storing the PNG files.
    width: float
        Width of plots in inches.
    height: float
        Height of plots in inches.
    dpi: int
        Dots per inch.
    f_start: float
        Starting frequency or None
    f_stop: float
        Stopping frequency or None
    """

    # Get directory path for storing PNG file
    if not isdir(plot_dir):
        os.mkdir(plot_dir)

    # Begin PNG file collection.
    png_collection = []

    # Generate a plot/PNG for each frequency chunk.
    logger.info("Image width = {}, height = {}, dpi = {}"
                .format(width, height, dpi))
    count = 0
    num_files = len(file_list)
    for input_file in file_list:

        # read in data
        wf = Waterfall(input_file, f_start=f_start, f_stop=f_stop)

        # Validate frequency range.
        count += 1
        logger.info("stax: Processing file {} of {}: {}"
                    .format(count, num_files, input_file))

        # Make plot for this frequency chunk with plot_waterfall().
        source_name = wf.header["source_name"]
        plt.figure(source_name, figsize=(width, height), dpi=dpi)
        plot_waterfall(wf, f_start=f_start, f_stop=f_stop)

        # Save the figures.
        path_png = plot_dir + "stax_" + str(count) + "_" + source_name + ".png"
        png_collection.append(path_png)
        plt.savefig(path_png, dpi=dpi)
        logger.info("Saved plot: {}".format(path_png))
        plt.close("all")

    return png_collection


def cmd_tool(args=None):
    r"""Coomand line parser"""
    parser = ArgumentParser(description="Make waterfall plots from a single file.")
    parser.add_argument('file_list', type=str, nargs='+', help='List of files to plot')
    parser.add_argument("--plot_dir", "-p", type=str, default=".", help="Directory to receive plots (.png).  Default: current directory.")
    parser.add_argument('--f_start', type=float, default=None, help='Start frequency.  Default: None.')
    parser.add_argument('--f_stop', type=float, default=None, help='Stop frequency.  Default: None.')
    parser.add_argument("--stitch", "-s", type=str, default="v", choices=["n", "v", "h"],
                        help="Stitch files: n(no stitching), v(vertical), or h(horizontal).  Default: v.")
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
    if args.dpi < 50:
        logger.error("stax: --dpi must be at least 50 but I saw: {}!".format(args.dpi))
        os.system("stax -h")
        sys.exit(86)
    if args.width < 6:
        logger.error("stax: --width must be at least 6 but I saw: {}!".format(args.width))
        os.system("stax -h")
        sys.exit(86)
    if args.height < 5:
        logger.error("stax: --height must be at least 5 but I saw: {}!".format(args.height))
        os.system("stax -h")
        sys.exit(86)

    # Make the plots.
    logger.info("Begin")
    args.plot_dir += "/"
    png_collection = make_waterfall_plots(args.file_list,
                                          args.plot_dir,
                                          args.width,
                                          args.height,
                                          args.dpi,
                                          args.f_start,
                                          args.f_stop)
    if args.stitch != "n":
        path_saved_png = args.plot_dir + "stax_stitched_" + args.stitch + ".png"
        image_stitch(args.stitch, png_collection, path_saved_png)
        for file in png_collection:
            os.remove(file)

    logger.info("End")


if __name__ == "__main__":
    cmd_tool()
