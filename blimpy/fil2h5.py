#!/usr/bin/env python
'''
    Simple script for quicly making an h5 file from a .fil.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)
'''

from blimpy import Waterfall
from optparse import OptionParser
import sys

def make_h5_file():
    ''' Converts file to h5 format. Saves output in current dir.
    '''

    p = OptionParser()
    p.set_usage('python fil2h5.py <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_option('-o', '--out_dir', dest='out_dir', type='str', default='./', help='Location for output files. Default: local dir. ')
    opts, args = p.parse_args(sys.argv[1:])

    if len(args)!=1:
        print 'Please specify a file name \nExiting.'
        sys.exit()
    else:
        filename = args[0]

    fil_file = Waterfall(filename)
    new_filename = opts.out_dir+filename.replace('.fil','.h5').split('/')[-1]
    fil_file.write_to_hdf5(new_filename)

if __name__ == "__main__":
    make_h5_file()
