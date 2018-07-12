"""
# test_voyager_data_load.py

The hard-coded numbers in these tests can be found in the voyager_test_setup.ipynb

"""

import blimpy as bl
import numpy as np
import pylab as plt


def test_max_data_array_size():
    fw = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.fil', max_data_array_size=1)
    fw = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.h5',  max_data_array_size=1)

if __name__ == "__main__":
    test_max_data_array_size()
