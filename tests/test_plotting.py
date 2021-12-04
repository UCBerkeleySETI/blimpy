import os
import blimpy as bl
import numpy as np
from pprint import pprint
import pylab as plt

from tests.data import voyager_fil, voyager_h5
from blimpy.plotting import plot_waterfall, plot_spectrum, plot_spectrum_min_max, \
    plot_kurtosis, plot_time_series, plot_all

TEST_DATA_DIR = os.path.dirname(voyager_h5)

def test_plot_waterfall():
    """ Load Voyager dataset and test plotting """

    a = bl.Waterfall(voyager_h5)

    plt.figure("TEST PLOTTING", figsize=(10, 8))
    plt.subplot(3, 2, 1)
    plot_waterfall(a)

    plt.subplot(3, 2, 2)
    plot_spectrum(a)

    plt.subplot(3, 2, 3)
    plot_spectrum_min_max(a)

    plt.subplot(3, 2, 4)
    plot_kurtosis(a)

    plt.subplot(3, 2, 5)
    plot_time_series(a)

    plt.tight_layout()
    plt.savefig(TEST_DATA_DIR + "/test_plotting.png")

    plt.figure("TEST PLOT_ALL", figsize=(10, 8))
    plot_all(a)
    plt.savefig(TEST_DATA_DIR + "/test_plotting_plot_all.png")


def test_plot_waterfall_classmethod():
    """ Load Voyager dataset and test plotting """

    a = bl.Waterfall(voyager_h5)

    plt.figure("TEST PLOTTING CLASS", figsize=(10, 8))
    plt.subplot(3, 2, 1)
    a.plot_waterfall()

    plt.subplot(3, 2, 2)
    a.plot_spectrum()

    plt.subplot(3, 2, 3)
    a.plot_spectrum_min_max()

    plt.subplot(3, 2, 4)
    a.plot_kurtosis()

    plt.subplot(3, 2, 5)
    a.plot_time_series()
    plt.tight_layout()

    plt.savefig(TEST_DATA_DIR + "/test_plotting_classmethod.png")

    plt.figure("TEST PLOT_ALL CLASS", figsize=(10, 8))
    a.plot_all()
    plt.savefig(TEST_DATA_DIR + "/test_plotting_plot_all_classmethod.png")


if __name__ == "__main__":
    test_plot_waterfall()
    test_plot_waterfall_classmethod()

