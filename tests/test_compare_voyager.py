import blimpy as bl
import numpy as np
from pprint import pprint

def compare_waterfall_fil_to_h5():
    """ Load Voyager dataset and test that both fil and hdf5 readers return same headers and data """

    print("Loading FIL and HDF5 data with Waterfall()..."),
    a = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.h5')
    b = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.fil')
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

if __name__ == "__main__":
    compare_waterfall_fil_to_h5()