#!/usr/bin/env python
"""
    Simple script for quickly making a .fil file from a .h5.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
"""

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


def bl_scrunch(filename, out_dir='./', new_filename=None, max_load=None, f_scrunch=None):
    """ Frequency scrunch (lower resolution by averaging) filename

    Args:
        filename (str): Name of file to open
        out_dir (str):
    """

    fil_file = Waterfall(filename, max_load=max_load)
    if not new_filename:
        new_filename = out_dir+filename.replace('.h5', '.scrunched.h5').split('/')[-1]

    print("Using fscrunch %i" % f_scrunch)
    fil_file.write_to_hdf5(new_filename, f_scrunch=f_scrunch)

def cmd_tool(args=None):
    """  Command line utility for converting HDF5 (.h5) to Sigproc filterbank (.fil) format

    Usage:
        h52fil <FULL_PATH_TO_FIL_FILE> [options]

    Options:
        -h, --help            show this help message and exit
        -o OUT_DIR, --out_dir=OUT_DIR
                              Location for output files. Default: local dir.
        -n NEW_FILENAME, --new_filename=NEW_FILENAME
                              New filename. Default: replaces extension to .fil
        -d, --delete_input    This option deletes the input file after conversion.
        -l MAX_LOAD           Maximum data limit to load. Default:1GB
    """

    p = OptionParser()
    p.set_usage('Command line utility for converting HDF5 (.h5) to Sigproc filterbank (.fil) format \n >>h52fil <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_option('-o', '--out_dir', dest='out_dir', type='str', default='./',
                 help='Location for output files. Default: local dir. ')
    p.add_option('-n', '--new_filename', dest='new_filename', type='str', default='',
                 help='New filename. Default: replaces extension to .scrunched.h5')
    p.add_option('-d', '--delete_input', dest='delete_input', action='store_true', default=False,
                 help='This option deletes the input file after conversion.')
    p.add_option('-f', '--fscrunch', dest='f_scrunch', default=1, type=int,
                 help='Average (aka scrunch) across frequency. Number of channels to average together.')
    p.add_option('-l', action='store', default=None, dest='max_load', type=float,
                 help='Maximum data limit to load. Default:1GB')

    if args is None:
        opts, args = p.parse_args(sys.argv[1:])
    else:
        opts, args = p.parse_args(args)

    if len(args) != 1:
        logger.info('Please specify a file name. \nExiting.')
        sys.exit()
    else:
        filename = args[0]

    if opts.f_scrunch == 1:
        logger.info('Please specify frequency scrunch amount with -f \nExiting.')
        sys.exit()

    bl_scrunch(filename, out_dir=opts.out_dir, new_filename=opts.new_filename,
               max_load=opts.max_load, f_scrunch=opts.f_scrunch)

    if opts.delete_input:
        logger.info("'Deleting input file: %s"%(filename))
        os.remove(filename)

if __name__ == "__main__":

    cmd_tool()
