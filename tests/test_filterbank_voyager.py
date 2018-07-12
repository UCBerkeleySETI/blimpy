import blimpy as bl
from blimpy import sigproc
import numpy as np
from pprint import pprint
import pytest

def compare_filterbank_fil_to_h5():
    """ Load Voyager dataset and test that both fil and hdf5 readers return same headers and data """

    print("Loading FIL and HDF5 data with Waterfall()..."),
    a = bl.Filterbank('Voyager_data/Voyager1.single_coarse.fine_res.h5')
    b = bl.Filterbank('Voyager_data/Voyager1.single_coarse.fine_res.fil')
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

def test_plotting_doesnt_cause_exceptions():
    """ Try running the plotting routines. They should not raise expections even without X windows """
    a = bl.Filterbank('Voyager_data/Voyager1.single_coarse.fine_res.h5')
    b = bl.Filterbank('Voyager_data/Voyager1.single_coarse.fine_res.fil')

    a.plot_all()
    a.plot_kurtosis()
    a.plot_spectrum()
    a.plot_spectrum_min_max()
    a.plot_waterfall()
    a.plot_time_series()

    b.plot_all()
    b.plot_kurtosis()
    b.plot_spectrum()
    b.plot_spectrum_min_max()
    b.plot_waterfall()
    b.plot_time_series()

def test_sigproc_is_fil():
    """ Check that the is_fil function works """

    assert sigproc.is_filterbank('Voyager_data/Voyager1.single_coarse.fine_res.h5') is False
    assert sigproc.is_filterbank('Voyager_data/Voyager1.single_coarse.fine_res.fil') is True

def test_file_wrapper_open_file():
    from blimpy.file_wrapper import open_file
    open_file('Voyager_data/Voyager1.single_coarse.fine_res.h5')
    open_file('Voyager_data/Voyager1.single_coarse.fine_res.fil')

    with pytest.raises(NotImplementedError):
        open_file('run_tests.sh')


if __name__ == "__main__":
    compare_filterbank_fil_to_h5()
    test_plotting_doesnt_cause_exceptions()
    test_sigproc_is_fil()
    test_file_wrapper_open_file()