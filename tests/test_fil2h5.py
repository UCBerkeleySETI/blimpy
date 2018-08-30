"""
# test_fil2h5
"""

import os

import blimpy as bl
from tests.data import voyager_fil


def test_fil2h5_conversion():
    '''Tests the conversion of fil files into h5 in both light and heavy modes.
    '''

    # Creating test file.
    bl.fil2h5.make_h5_file(voyager_fil, new_filename='test.h5')

    # Creating a "large" test file.
    bl.fil2h5.make_h5_file(voyager_fil, new_filename='test_large.h5', max_load=0.001)

    # Testing filename
    bl.fil2h5.make_h5_file(voyager_fil, new_filename='test')

    # Deleting test file
    os.remove('test.h5')
    os.remove('test_large.h5')

if __name__ == "__main__":
    test_fil2h5_conversion()
