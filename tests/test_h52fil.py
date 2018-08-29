"""
# test_h52fil
"""

import blimpy as bl
import numpy as np
import pylab as plt
import os


def test_h52fil_conversion():
    ''' Tests the conversion of fil files into h5 in both light and heavy modes.
    '''

    #Creating test file.
    bl.h52fil.make_fil_file('Voyager_data/Voyager1.single_coarse.fine_res.h5', new_filename = 'test.fil')

    #Creating a "large" test file.
    bl.h52fil.make_fil_file('Voyager_data/Voyager1.single_coarse.fine_res.h5', new_filename = 'test_large.fil', max_load = 0.001)

    #Deleting test file
    os.remove('test.fil')
    os.remove('test_large.fil')

if __name__ == "__main__":
    test_h52fil_conversion()


