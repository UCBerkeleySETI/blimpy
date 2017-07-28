#!/usr/bin/env python
'''
    Simple script for quicly making an h5 file from a .fil.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
'''

from .waterfall import Waterfall
from optparse import OptionParser
import sys

def make_h5_file(filename,out_dir='./'):
    ''' Converts file to HDF5 (.h5) format. Default saves output in current dir.
    '''

    fil_file = Waterfall(filename)
    new_filename = out_dir+filename.replace('.fil','.h5').split('/')[-1]
    fil_file.write_to_hdf5(new_filename)

if __name__ == "__main__":

    p = OptionParser()
    p.set_usage('Command line utility for converting Sigproc filterbank (.fil) to  HDF5 (.h5) format  \n >>fil2h5 <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_option('-o', '--out_dir', dest='out_dir', type='str', default='./', help='Location for output files. Default: local dir. ')
    opts, args = p.parse_args(sys.argv[1:])

    if len(args)!=1:
        print('Please specify a file name \nExiting.')
        sys.exit()
    else:
        filename = args[0]

    make_h5_file(filename, out_dir = opts.out_dir)
