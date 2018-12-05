import blimpy as bl
import numpy as np
from pprint import pprint
from tests.data import voyager_fil, voyager_h5
from blimpy.plotting.config import plt

def compare_waterfall_fil_to_h5():
    """ Load Voyager dataset and test that both fil and hdf5 readers return same headers and data """

    print("Loading FIL and HDF5 data with Waterfall()..."),
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)
    print("OK")

    print("Reading headers..")
    print("\nHDF5 file header:")
    pprint(a.header)
    print("\nFIL file header:")
    pprint(b.header)
    print("Headers are loading OK")

    print("\nChecking header values match..."),
    for key in b.header.keys():
        assert b.header[key] == a.header[key]
    print("OK")

    print("Checking datatype matches..."),
    assert a.data.dtype == b.data.dtype
    print("OK")

    print("Checking data matches..."),
    assert np.allclose(a.data, b.data)
    assert a.data.dtype == b.data.dtype
    print("OK")

def compare_waterfall_fil_to_h5_methods_and_attributes():
    """ Compare attributes and check methods """
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)

    print("Comparing attributes of classes match where expected")
    assert a.beam_axis == b.beam_axis
    assert a.freq_axis == b.freq_axis
    assert a.time_axis == b.time_axis
    #assert a.stokes_axis == b.stokes_axis  # Fil shouldn't have stokes axis ...

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

    # TODO: 'compute_lsrk' -- need PYSLALIB
    # TODO: 'compute_lst' -- need PYSLALIB

    # Expected to differ between: 'ext', 'file_size_bytes', 'filename', 'container'
    # Unused? Inherited from Filterbank.py 'gen_from_header', 'generate_freqs', 'read_filterbank', 'read_hdf5'

    # TODO: 'grab_data',
    # TODO: 'read_data',
    # TODO: 'write_to_fil',
    # TODO: 'write_to_filterbank',
    # TODO: 'write_to_hdf5']

    dir_a = dir(a)
    dir_b = dir(b)

    print("Attr/methods in HDF5 but not in FIL:")
    for item in dir_a:
        if item not in dir_b:
            print(item)

    print("Attr/methods in FIL but not in HDF5:")
    for item in dir_b:
        if item not in dir_a:
            print(item)

def compare_waterfall_fil_h5_conatiners():
    """ Compare the two containers for fil.container and h5.container """
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)

    dir_a = dir(a.container)
    dir_b = dir(b.container)

    print("Attr/methods in HDF5 container but not in FIL:")
    for item in dir_a:
        if item not in dir_b:
            print(item)

    print("Attr/methods in FIL container but not in HDF5:")
    for item in dir_b:
        if item not in dir_a:
            print(item)

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

if __name__ == "__main__":
    compare_waterfall_fil_to_h5()
    compare_waterfall_fil_to_h5_methods_and_attributes()
    compare_waterfall_fil_h5_conatiners()
    test_plotting_doesnt_cause_exceptions()