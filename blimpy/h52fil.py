#!/usr/bin/env python
"""
    Simple script for making a .fil file from a .h5.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
"""


from argparse import ArgumentParser
import sys
import os
from blimpy import Waterfall
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

    parser = ArgumentParser('Command line utility for converting HDF5 (.h5) to Sigproc filterbank (.fil) format \n >>h52fil <FULL_PATH_TO_FIL_FILE> [options]')
    parser.add_argument("filepath_in", type=str, help="Path of input HDF5 Filterbank file")
    parser.add_argument('-o', '--out_dir', dest='out_dir', type=str, default='./',
                 help='Location for output files. Default: local dir. ')
    parser.add_argument('-n', '--new_filename', dest='new_filename', type=str, default='',
                 help='New filename. Default: replaces extension to .fil')
    parser.add_argument('-d', '--delete_input', dest='delete_input', action='store_true', default=False,
                 help='This option deletes the input file after conversion.')
    parser.add_argument('-l', action='store', default=None, dest='max_load', type=float,
                 help='Maximum data limit to load. Default:1GB')
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    make_fil_file(args.filepath_in, out_dir = args.out_dir, new_filename=args.new_filename, max_load = args.max_load)

    if args.delete_input:
        logger.info("Deleting input file: {}".format(args.filepath_in))
        os.remove(args.filepath_in)

if __name__ == "__main__":

    cmd_tool()
