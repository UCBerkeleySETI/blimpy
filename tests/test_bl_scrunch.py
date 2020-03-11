"""
# test_h52fil
"""

import os

import blimpy as bl
from tests.data import voyager_h5


def test_scrunch():
    ''' Tests the conversion of fil files into h5 in both light and heavy modes.
    '''

    # Creating test file.
    bl.bl_scrunch.bl_scrunch(voyager_h5, new_filename='test.scrunched.h5', f_scrunch=8)

    # Deleting test file
    os.remove('test.scrunched.h5')

if __name__ == "__main__":
    test_scrunch()
