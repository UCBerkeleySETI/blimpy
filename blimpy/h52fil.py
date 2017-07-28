#!/usr/bin/env python
'''
    Simple script for quicly making an .fil file from a .h5.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
'''

from .waterfall import Waterfall
from optparse import OptionParser
import sys

def make_fil_file(filename,out_dir='./'):
    ''' Converts file to Sigproc filterbank (.fil) format.  Default saves output in current dir.
    '''

    fil_file = Waterfall(filename)
    new_filename = out_dir+filename.replace('.h5','.fil').split('/')[-1]
    fil_file.write_to_fil(new_filename)

if __name__ == "__main__":

    p = OptionParser()
    p.set_usage('Command line utility for converting HDF5 (.h5) to Sigproc filterbank (.fil) format \n >>h52fil <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_option('-o', '--out_dir', dest='out_dir', type='str', default='./', help='Location for output files. Default: local dir. ')
    opts, args = p.parse_args(sys.argv[1:])

    if len(args)!=1:
        print('Please specify a file name \nExiting.')
        sys.exit()
    else:
        filename = args[0]

    make_fil_file(filename, out_dir = opts.out_dir)
