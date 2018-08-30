#!/usr/bin/env python
'''
    Simple script for quicly making an h5 file from a .fil.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
'''

try:
    from .waterfall import Waterfall
except:
    from waterfall import Waterfall

from optparse import OptionParser
import sys
import os

#------
# Logging set up
import logging
logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    format = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    format = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'

logging.basicConfig(format=format,stream=stream,level = level_log)
#------


def make_h5_file(filename,out_dir='./', new_filename = None, max_load = None):
    ''' Converts file to HDF5 (.h5) format. Default saves output in current dir.
    '''

    fil_file = Waterfall(filename, max_load = max_load)
    if not new_filename:
        new_filename = out_dir+filename.replace('.fil','.h5').split('/')[-1]

    if '.h5' not in new_filename:
        new_filename = new_filename+'.h5'

    fil_file.write_to_hdf5(new_filename)

def cmd_tool():
    '''
    '''

    p = OptionParser()
    p.set_usage('Command line utility for converting Sigproc filterbank (.fil) to  HDF5 (.h5) format  \n >>fil2h5 <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_option('-o', '--out_dir', dest='out_dir', type='str', default='./', help='Location for output files. Default: local dir. ')
    p.add_option('-n', '--new_filename', dest='new_filename', type='str', default='', help='New name. Default: replaces extention to .h5')
    p.add_option('-d', '--delete_input', dest='delete_input', action='store_true', default=False, help='This option deletes the input file after conversion.')
    p.add_option('-l', action='store', default=None, dest='max_load', type=float,help='Maximum data limit to load. Default:1GB')

    opts, args = p.parse_args(sys.argv[1:])

    if len(args)!=1:
        logger.info('Please specify a file name \nExiting.')
        sys.exit()
    else:
        filename = args[0]

    make_h5_file(filename, out_dir = opts.out_dir, new_filename = opts.new_filename, max_load = opts.max_load)

    if opts.delete_input:
        logger.info("'Deleting input file: %s"%(filename))
        os.remove(filename)

if __name__ == "__main__":

    cmd_tool()
