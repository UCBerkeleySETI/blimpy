"""
Very small module with one test: test_write_to_fil()
"""
import os
import blimpy as bl

from tests.data import voyager_h5

OUTDIR = os.path.dirname(voyager_h5) + "/"

def test_write_to_fil():
    """ Load Voyager dataset and test plotting """

    a = bl.Waterfall(voyager_h5)
    a.write_to_fil(OUTDIR + 'test_out.fil')

if __name__ == "__main__":
    test_write_to_fil()
