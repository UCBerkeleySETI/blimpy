import blimpy as bl
import numpy as np
from pprint import pprint
import pytest
from tests.data import voyager_fil, voyager_h5
from blimpy.plotting.config import plt


def compare_filterbank_fil_to_h5():
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


def test_plotting_doesnt_cause_exceptions():
    """ Try running the plotting routines. They should not raise expections even without X windows """
    a = bl.Waterfall(voyager_h5)
    b = bl.Waterfall(voyager_fil)

    a.plot_all()
    a.plot_kurtosis()
    a.plot_spectrum()
    a.plot_spectrum_min_max()
    a.plot_waterfall()
    a.plot_time_series()
    plt.clf() # Fix issue #140

    b.plot_all()
    b.plot_kurtosis()
    b.plot_spectrum()
    b.plot_spectrum_min_max()
    b.plot_waterfall()
    b.plot_time_series()
    plt.clf()


def test_cmdtool():
    with pytest.raises(SystemExit):
        bl.waterfall.cmd_tool(args=[])


if __name__ == "__main__":
    compare_filterbank_fil_to_h5()
    test_plotting_doesnt_cause_exceptions()

    test_cmdtool()
