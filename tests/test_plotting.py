import blimpy as bl
import numpy as np
from pprint import pprint
import pylab as plt

from data import voyager_fil, voyager_h5
from blimpy.plotting import plot_waterfall, plot_spectrum, plot_spectrum_min_max, \
    plot_kurtosis, plot_time_series, plot_all


def test_plot_waterfall():
    """ Load Voyager dataset and test plotting """

    print("Loading FIL and HDF5 data with Waterfall()..."),
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)
    print("OK")

    plt.figure(figsize=(10, 8))
    plt.subplot(3,2,1)
    plot_waterfall(a)

    plt.subplot(3,2,2)
    plot_spectrum(a)

    plt.subplot(3,2,3)
    plot_spectrum_min_max(a)

    plt.subplot(3,2,4)
    plot_kurtosis(a)

    plt.subplot(3,2,5)
    plot_time_series(a)
    plt.tight_layout()
    plt.show()

    plot_all(a)
    plt.show()



if __name__ == "__main__":
    test_plot_waterfall()