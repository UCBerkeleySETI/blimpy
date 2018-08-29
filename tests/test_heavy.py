"""
# test_heavy.py

"""

import blimpy as bl
from tests.data import voyager_fil, voyager_h5


def test_max_data_array_size():
    fw = bl.Waterfall(voyager_fil, max_load=0.001)
    fw = bl.Waterfall(voyager_h5,  max_load=0.001)

if __name__ == "__main__":
    test_max_data_array_size()
