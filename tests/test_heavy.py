"""
# test_heavy.py

"""

import blimpy as bl
import numpy as np
import pylab as plt


def test_max_data_array_size():
    fw = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.fil', max_load=0.001)
    fw = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.h5',  max_load=0.001)

if __name__ == "__main__":
    test_max_data_array_size()
