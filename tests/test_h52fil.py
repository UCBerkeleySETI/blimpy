"""
# test_h52fil
"""

import os

import blimpy as bl
from tests.data import voyager_h5


def test_h52fil_conversion():
    ''' Tests the conversion of fil files into h5 in both light and heavy modes.
    '''

    # Creating test file.
    bl.h52fil.make_fil_file(voyager_h5, new_filename='test.fil')

    # Creating a "large" test file.
    bl.h52fil.make_fil_file(voyager_h5, new_filename='test_large.fil', max_load=0.001)

    # Testing filename
    bl.h52fil.make_fil_file(voyager_h5, new_filename='test')

    # Deleting test file
    os.remove('test.fil')
    os.remove('test_large.fil')

if __name__ == "__main__":
    test_h52fil_conversion()
