from blimpy import sigproc
import blimpy as bl
from tests.data import voyager_fil, voyager_h5
import numpy as np
import os


def test_sigproc_is_fil():
    """ Check that the is_fil function works """

    assert sigproc.is_filterbank(voyager_h5) is False
    assert sigproc.is_filterbank(voyager_fil) is True


def test_sigproc_generate_headers():
    """ Test if you can generate headers OK from files """
    a = bl.Filterbank(voyager_h5)
    b = bl.Filterbank(voyager_fil)
    sigproc.generate_sigproc_header(a)
    sigproc.generate_sigproc_header(b)

def test_fil_write():
    try:
        a = bl.Filterbank(voyager_h5)
        b = bl.Filterbank(voyager_fil)

        a.write_to_filterbank('test.fil')
        b.write_to_filterbank('test2.fil')

        c = bl.Filterbank('test.fil')
        d = bl.Filterbank('test2.fil')

        for key in a.header.keys():
            if key != b'DIMENSION_LABELS':
                assert a.header[key] == c.header[key]
                assert key in c.header.keys()
                assert a.header[key] == d.header[key]
                assert key in d.header.keys()

        assert np.allclose(a.data, c.data)
        assert np.allclose(a.data, d.data)
    except AssertionError:
        print(key, a.header[key], b.header[key], c.header[key], d.header[key])
        raise

    finally:
        os.remove('test.fil')
        os.remove('test2.fil')

if __name__ == "__main__":
    test_sigproc_is_fil()
    test_sigproc_generate_headers()
    test_fil_write()
