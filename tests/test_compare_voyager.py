import blimpy as bl
import numpy as np
from pprint import pprint
from tests.data import voyager_fil, voyager_h5
from blimpy.plotting.config import plt


def test_compare_waterfall_fil_to_h5():
    """ Load Voyager dataset and test that both fil and hdf5 readers return same headers and data """

    print("Loading FIL and HDF5 data with Waterfall()..."),
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)

    print("Reading headers..")
    pprint(a.header)
    print("\nFIL file header:")
    pprint(b.header)
    print("Headers are loading OK")

    print("\nChecking header values match..."),
    for key in b.header.keys():
        assert b.header[key] == a.header[key]

    print("Checking datatype matches..."),
    assert a.data.dtype == b.data.dtype

    print("Checking data matches..."),
    assert np.allclose(a.data, b.data)
    assert a.data.dtype == b.data.dtype


def test_waterfall_fil_to_h5_methods_and_attributes():
    """ Compare attributes and check methods """
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)

    print("Comparing attributes of classes match where expected")
    assert a.beam_axis == b.beam_axis
    assert a.freq_axis == b.freq_axis
    assert a.time_axis == b.time_axis

    assert a.calc_n_coarse_chan() == b.calc_n_coarse_chan()
    assert a.file_shape == b.file_shape

    assert a.n_channels_in_file == b.n_channels_in_file
    assert a.n_ints_in_file == b.n_ints_in_file
    assert a.selection_shape == b.selection_shape

    print("Checking if basic methods run without raising Exceptions")
    # Check they can be run
    a.container.populate_freqs()
    a.container.populate_timestamps()
    a.info()
    a.blank_dc(1)
    a.calibrate_band_pass_N1()

    b.container.populate_freqs()
    b.container.populate_timestamps()
    b.info()
    b.blank_dc(1)
    b.calibrate_band_pass_N1()

    dir_a = dir(a)
    dir_b = dir(b)

    print("Attr/methods in HDF5 but not in FIL:")
    for item in dir_a:
        if item not in dir_b:
            raise ValueError("HDF5 item is not in FIL:" + str(item))

    print("Attr/methods in FIL but not in HDF5:")
    for item in dir_b:
        if item not in dir_a:
            raise ValueError("FIL item is not in HDF5:" + str(item))


def test_plotting_doesnt_cause_exceptions():
    """ Try running the plotting routines. They should not raise expections even without X windows """
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)

    a.plot_all()
    plt.clf()
    a.plot_kurtosis()
    plt.clf()
    a.plot_spectrum()
    plt.clf()
    a.plot_spectrum_min_max()
    plt.clf()
    a.plot_waterfall()
    plt.clf()
    a.plot_time_series()
    plt.clf()

    b.plot_all()
    plt.clf()
    b.plot_kurtosis()
    plt.clf()
    b.plot_spectrum()
    plt.clf()
    b.plot_spectrum_min_max()
    plt.clf()
    b.plot_waterfall()
    plt.clf()
    b.plot_time_series()
    plt.clf()

