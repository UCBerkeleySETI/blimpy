#!/usr/bin/env python3

r"""
Make waterfall plots of a file set, view from top to bottom.
"""

import os
import sys
from os.path import dirname, abspath, isdir
import gc
from argparse import ArgumentParser
import logging
logger_name = "stax"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.INFO)

# Plotting packages import
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("agg")

# Math/Science package imports
import numpy as np

# Blimpy imports
import blimpy as bl
from blimpy.utils import rebin

# Plotting constants
fontsize=16
font = {"family" : "DejaVu Sans",
"size" : fontsize}
MAX_IMSHOW_POINTS = (4096, 1268)


def plot_waterfall(wf, f_start=None, f_stop=None, **kwargs):
    r"""
    Plot waterfall of data in a .fil or .h5 file.

    Parameters
    ----------
    wf : blimpy.Waterfall object
        Waterfall object of an H5 or Filterbank file containing the dynamic spectrum data.
    f_start : float
        Start frequency, in MHz.
    f_stop : float
        Stop frequency, in MHz.
    kwargs : dict
        Keyword args to be passed to matplotlib imshow().

    Notes
    -----
    Plot a single-panel waterfall plot (frequency vs. time vs. intensity)
    for one of the files in the set of interest, at the
    frequency of the expected event.
    """

    # Extract source name.
    source_name = wf.header["source_name"]

    # prepare font
    matplotlib.rc("font", **font)

    # Load in the data from fil
    plot_f, plot_data = wf.grab_data(f_start=f_start, f_stop=f_stop)

    # Make sure waterfall plot is under 4k*4k
    dec_fac_x, dec_fac_y = 1, 1

    # rebinning data to plot correctly with fewer points
    try:
        if plot_data.shape[0] > MAX_IMSHOW_POINTS[0]:
            dec_fac_x = plot_data.shape[0] / MAX_IMSHOW_POINTS[0]
        if plot_data.shape[1] > MAX_IMSHOW_POINTS[1]:
            dec_fac_y =  int(np.ceil(plot_data.shape[1] /  MAX_IMSHOW_POINTS[1]))
        plot_data = rebin(plot_data, dec_fac_x, dec_fac_y)
    except Exception as ex:
        print("\n*** Oops, grab_data returned plot_data.shape={}, plot_f.shape={}"
              .format(plot_data.shape, plot_f.shape))
        print("Waterfall info for {}:".format(wf.filename))
        wf.info()
        raise ValueError("*** Something is wrong with the grab_data output!") from ex

    # Rolled back PR #82

    # determine extent of the plotting panel for imshow
    extent=(plot_f[0], plot_f[-1], (wf.timestamps[-1]-wf.timestamps[0])*24.*60.*60, 0.0)

    # plot and scale intensity (log vs. linear)
    kwargs["cmap"] = kwargs.get("cmap", "viridis")
    plot_data = 10.0 * np.log10(plot_data)

    # get normalization parameters
    vmin = plot_data.min()
    vmax = plot_data.max()
    normalized_plot_data = (plot_data - vmin) / (vmax - vmin)

    # display the waterfall plot
    this_plot = plt.imshow(normalized_plot_data,
        aspect="auto",
        rasterized=True,
        interpolation="nearest",
        extent=extent,
        **kwargs
    )

    # add plot labels
    plt.xlabel("Frequency [Hz]",fontdict=font)
    plt.ylabel("Time [s]",fontdict=font)

    # add source name
    ax = plt.gca()
    plt.text(0.03, 0.8, source_name, transform=ax.transAxes, bbox=dict(facecolor="white"))

    return this_plot


def sort2(x, y):
    r""" Return lowest value, highest value"""
    if y < x:
        return y, x
    return x, y


def make_waterfall_plots(file_list, plot_dir, f_start=None, f_stop=None, **kwargs):
    r"""
    Make waterfall plots of a file set, view from top to bottom.

    Parameters
    ----------
    file_list : str
        List of filterbank files to plot in a stacked mode.
    f_start : float
        Start frequency, in MHz.
    f_stop : float
        Stop frequency, in MHz.
    source_name_list : list
        List of source names in the set, in order.
    kwargs : dict
        Keyword args to be passed to matplotlib imshow().
    """

    # prepare for plotting
    matplotlib.rc("font", **font)

    # set up the sub-plots
    n_plots = len(file_list)
    fig = plt.subplots(n_plots, sharex=True, sharey=True,figsize=(10, 2*n_plots))

    # get directory path for storing PNG files
    if plot_dir is None:
        dirpath = dirname(abspath(file_list[0])) + "/"
    else:
        if not isdir(plot_dir):
            os.mkdir(plot_dir)
        dirpath = plot_dir


    # read in data for the first panel
    max_load = bl.calcload.calc_max_load(file_list[0])
    #print("plot_event make_waterfall_plots: max_load={} is required for {}".format(max_load, file_list[0]))
    wf1 = bl.Waterfall(file_list[0], max_load=max_load)
    t0 = wf1.header["tstart"]
    source_name = wf1.header["source_name"]

    # Fix frequency boundaries if required.
    freqs = wf1.container.populate_freqs()
    if f_start is None:
        f_start = freqs[0]
    if f_stop is None:
        f_stop = freqs[-1]
    the_lowest, the_highest = sort2(f_start, f_stop)

    # rebin data to plot correctly with fewer points
    plot_f1, plot_data1 = wf1.grab_data(f_start=the_lowest, f_stop=the_highest)
    dec_fac_x, dec_fac_y = 1, 1
    if plot_data1.shape[0] > MAX_IMSHOW_POINTS[0]:
        dec_fac_x = plot_data1.shape[0] / MAX_IMSHOW_POINTS[0]
    if plot_data1.shape[1] > MAX_IMSHOW_POINTS[1]:
        dec_fac_y =  int(np.ceil(plot_data1.shape[1] /  MAX_IMSHOW_POINTS[1]))
    plot_data1 = rebin(plot_data1, dec_fac_x, dec_fac_y)

    # Compute the midpoint.
    the_midpoint = np.abs(the_lowest + the_highest) / 2.

    # Protect RAM utilisation.
    del wf1, freqs, plot_f1, plot_data1
    gc.collect()

    # Fill in each subplot for the full plot
    subplots = []
    for ii, filename in enumerate(file_list):
        logger.debug("make_waterfall_plots: file {} in list: {}".format(ii, filename))
        # identify panel
        subplot = plt.subplot(n_plots, 1, ii + 1)
        subplots.append(subplot)

        # read in data
        max_load = bl.calcload.calc_max_load(filename)
        wf = bl.Waterfall(filename, max_load=max_load)

        # Validate frequency range.
        freqs = wf.container.populate_freqs()
        ii_lowest, ii_highest = sort2(freqs[0], freqs[-1])
        logger.info("Processing: {}, freq lowest={}, highest={}".format(filename, ii_lowest, ii_highest))
        if the_lowest < ii_lowest or the_highest > ii_highest:
            logger.warning("Frequency range not compatible!  Ignoring this file.")
            # Protect RAM utilisation.
            del wf, freqs
            gc.collect()
            continue # Skip this file.

        # make plot with plot_waterfall
        source_name = wf.header["source_name"]
        this_plot = plot_waterfall(wf,
                                   f_start=the_lowest,
                                   f_stop=the_highest,
                                   **kwargs)

        # Title the full plot if processing the first file.
        if ii == 0:
            plot_title = "%s \n MJD:%5.5f (first file)" % (source_name, t0)
            plt.title(plot_title)

        # Format full plot.
        if ii < len(file_list)-1:
            plt.xticks(np.linspace(the_lowest, the_highest, num=4), ["","","",""])

        # Protect RAM utilisation.
        del wf, freqs
        gc.collect()

    # More overall plot formatting, axis labelling.
    factor = 1e6
    units = "Hz"
    xloc = np.linspace(the_lowest, the_highest, 5)
    xticks = [round(loc_freq) for loc_freq in (xloc - the_midpoint) * factor]
    if np.max(xticks) > 1000:
        xticks = [xt / 1000 for xt in xticks]
        units = "kHz"
    plt.xticks(xloc, xticks)
    plt.xlabel("Relative Frequency [%s] from %f MHz" % (units, the_midpoint), fontdict=font)

    # Add colorbar.
    cax = fig[0].add_axes([0.94, 0.11, 0.03, 0.77])
    fig[0].colorbar(this_plot, cax=cax, label="Normalized Power (Arbitrary Units)")

    # Adjust plots
    plt.subplots_adjust(hspace=0,wspace=0)

    # Save the figures.
    path_png = dirpath + source_name + "_fstart_{:0.6f}".format(f_start) \
                + "_fstop_{:0.6f}".format(f_stop) + ".png"
    plt.savefig(path_png, bbox_inches="tight")
    logger.info("Saved plot: {}".format(path_png))

    # show figure before closing if this is an interactive context
    mplbe = matplotlib.get_backend()
    logger.debug("make_waterfall_plots: backend = {}".format(mplbe))
    if mplbe != "agg":
        plt.show()

    # close all figure windows
    plt.close("all")


def cmd_line(args=None):
    r"""Coomand line parser"""
    parser = ArgumentParser(description='Make waterfall plots from a set of files, viewed from top to bottom.')
    parser.add_argument('file_list', type=str, nargs='+', help='List of files to plot')
    parser.add_argument('-f_start', type=float, default=None, help='Start frequency.  Default: None.')
    parser.add_argument('-f_stop', type=float, default=None, help='Stop frequency.  Default: None.')
    parser.add_argument('-plot_dir', type=str, default=".", help='Directory to receive the plot (.png).  Default: current directory.')
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    if args.file_list == []:
        os.system("stax -h")
        sys.exit(0)

    # Make the plots.
    args.plot_dir += "/"
    args.file_list = sorted(args.file_list)
    logger.info("Begin")
    make_waterfall_plots(args.file_list, args.plot_dir, args.f_start, args.f_stop)
    logger.info("End")


if __name__ == "__main__":
    cmd_line()
