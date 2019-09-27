import blimpy as bl
import numpy as np
from pprint import pprint
import pylab as plt

from tests.data import voyager_fil, voyager_h5

from blimpy.ephemeris import compute_lst, compute_lsrk

def test_compute_lst():
    """ Load Voyager dataset and test plotting """
    print("Loading HDF5 data with Waterfall()..."),
    a = bl.Waterfall(voyager_h5)
    print(compute_lst(a))


def test_compute_lsrk():
    a = bl.Waterfall(voyager_h5)
    print(compute_lsrk(a))


if __name__ == "__main__":
    test_compute_lst()
    test_compute_lsrk()