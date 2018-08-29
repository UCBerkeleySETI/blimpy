"""
# test_fil2h5
"""

import blimpy as bl
import numpy as np
import pylab as plt
import os


def test_fil2h5_conversion():
    ''' Tests the conversion of fil files into h5 in both light and heavy modes.
    '''

    #Creating test file.
    bl.fil2h5.make_h5_file('Voyager_data/Voyager1.single_coarse.fine_res.fil', new_filename = 'test.h5')

    #Creating a "large" test file.
    bl.fil2h5.make_h5_file('Voyager_data/Voyager1.single_coarse.fine_res.fil', new_filename = 'test_large.h5', max_load = 0.001)

    #Deleting test file
    os.remove('test.h5')
    os.remove('test_large.h5')

if __name__ == "__main__":
    test_fil2h5_conversion()


