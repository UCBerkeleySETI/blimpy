#!/usr/bin/env python
"""
    Simple script for making a .fil file from a .h5.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
"""

try:
    from .waterfall import Waterfall
except:
    from waterfall import Waterfall

from argparse import ArgumentParser
import sys
import os
from .utils import change_the_ext

import logging
logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    fmt = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    fmt = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'

logging.basicConfig(format=fmt, stream=stream, level = level_log)


def make_fil_file(filename,out_dir='./', new_filename=None, max_load = None):
    """ Converts file to Sigproc filterbank (.fil) format.  Default saves output in current dir.
    """

    wf = Waterfall(filename, max_load = max_load)

    if not new_filename:
        new_filename = out_dir + change_the_ext(filename, 'h5', 'fil').split('/')[-1]

    wf.write_to_fil(new_filename)


def cmd_tool(flags=None):
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

    p = ArgumentParser('Command line utility for converting HDF5 (.h5) to Sigproc filterbank (.fil) format \n >>h52fil <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_argument('-o', '--out_dir', dest='out_dir', type=str, default='./',
                 help='Location for output files. Default: local dir. ')
    p.add_argument('-n', '--new_filename', dest='new_filename', type=str, default='',
                 help='New filename. Default: replaces extension to .fil')
    p.add_argument('-d', '--delete_input', dest='delete_input', action='store_true', default=False,
                 help='This option deletes the input file after conversion.')
    p.add_argument('-l', action='store', default=None, dest='max_load', type=float,
                 help='Maximum data limit to load. Default:1GB')
    if flags is None:
        opts, args = p.parse_args(sys.argv[1:])
    else:
        opts, args = p.parse_args(flags)

    if len(args) != 1:
        logger.info('Please specify a file name \nExiting.')
        sys.exit()
    else:
        filename = args[0]

    make_fil_file(filename, out_dir = opts.out_dir, new_filename=opts.new_filename, max_load = opts.max_load)

    if opts.delete_input:
        logger.info("'Deleting input file: %s"%(filename))
        os.remove(filename)

if __name__ == "__main__":

    cmd_tool()
