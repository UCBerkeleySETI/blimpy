import blimpy as bl

from tests.data import voyager_h5

def test_write_to_fil():
    """ Load Voyager dataset and test plotting """

    a = bl.Waterfall(voyager_h5)
    a.write_to_fil('test_out.fil')

if __name__ == "__main__":
    test_write_to_fil()