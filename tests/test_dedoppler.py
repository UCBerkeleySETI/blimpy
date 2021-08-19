from os.path import dirname
import numpy as np
from astropy import units as u
import setigen as stg
import matplotlib.pyplot as plt
from blimpy.signal_processing.dedoppler import dedoppler_1
from blimpy import Waterfall
from blimpy.plotting import plot_waterfall
from tests.data import voyager_fil


PLOT_DIR = dirname(voyager_fil)
FIL_FILE = PLOT_DIR + "/test_dedoppler.fil"
PNG_FILE = PLOT_DIR + "/test_dedoppler.png"


# Plotting constants
fontsize = 16
font_dict = {"family" : "DejaVu Sans", "size" : fontsize}
N_PLOTS = 6


def sort2(x, y):
    r""" Return lowest value, highest value"""
    if y < x:
        return y, x
    return x, y


def plotter(counter, drift_rate):
    wf = Waterfall(FIL_FILE)
    dedoppler_1(wf, drift_rate)
    wf.header["source_name"] = "Dedoppler at D.R. {} Hz".format(drift_rate)
    plt.subplot(N_PLOTS, 1, counter)
    return plot_waterfall(wf)


def test_dedoppler_1():

    # Generate the Filterbank file.
    print("test_dedoppler_1: Creating Filterbank file {}".format(FIL_FILE))
    frame = stg.Frame(fchans=1024*u.pixel,
                      tchans=32*u.pixel,
                      df=2.7939677238464355*u.Hz,
                      dt=18.253611008*u.s,
                      fch1=6095.214842353016*u.MHz,
                      ascending=True)
    frame.add_noise(x_mean=10, noise_type='chi2')
    frame.add_signal(stg.constant_path(f_start=frame.get_frequency(index=200),
                                       drift_rate=2*u.Hz/u.s),
                                       stg.constant_t_profile(level=frame.get_intensity(snr=30)),
                                       stg.gaussian_f_profile(width=40*u.Hz),
                                       stg.constant_bp_profile(level=1))
    frame.save_fil(FIL_FILE)

    # Load Filterban file.
    print("test_dedoppler_1: Loading Filterbank file {}".format(FIL_FILE))
    wf = Waterfall(FIL_FILE)
    freqs = wf.get_freqs()
    the_lowest, the_highest = sort2(freqs[0], freqs[-1])
    the_midpoint = np.abs(the_lowest + the_highest) / 2.

    # Initialise plotting.
    print("test_dedoppler_1: Plotting to file {}".format(PNG_FILE))
    plt.subplots(N_PLOTS, sharex=True, sharey=True, figsize=(10, 2 * N_PLOTS))
    wf.header["source_name"] = "Initial Data"

    # Plot 1.
    plt.subplot(N_PLOTS, 1, 1)
    plot_waterfall(wf)

    # Plot #2.
    plotter(2, 1.0)

    # Plot #3.
    plotter(3, 1.5)

    # Plot #4.
    plotter(4, 2.2)

    # Plot #5.
    plotter(5, 2.7)

    # Plot #6.
    plotter(6, 3.6)

    # Finish up plots.
    plt.xticks(np.linspace(the_lowest, the_highest, num=4), ["","","",""])
    factor = 1e6
    units = "Hz"
    xloc = np.linspace(the_lowest, the_highest, 5)
    xticks = [round(loc_freq) for loc_freq in (xloc - the_midpoint) * factor]
    if np.max(xticks) > 1000:
        xticks = [xt / 1000 for xt in xticks]
        units = "kHz"
    plt.xticks(xloc, xticks)
    plt.xlabel("Relative Frequency [%s] from %f MHz" % (units, the_midpoint), fontdict=font_dict)
    plt.subplots_adjust(hspace=0, wspace=0)

    plt.savefig(PNG_FILE, dpi=200, bbox_inches="tight")
    print("test_dedoppler_1: End")


if __name__ == "__main__":
    test_dedoppler_1()
