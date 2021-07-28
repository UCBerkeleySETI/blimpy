#!/usr/bin/env python
"""
    Simple script for quicly making an h5 file from a .fil.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
"""

import sys
import os
from argparse import ArgumentParser

# Logging set up
import logging
logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    fmt = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    fmt = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'
logging.basicConfig(format=fmt,stream=stream,level = level_log)


from blimpy import Waterfall
from blimpy.io.hdf_writer import __write_to_hdf5_heavy as write_to_h5


def make_h5_file(filename, out_dir='./', new_filename=None, t_start=None, t_stop=None):
    """ Converts file to HDF5 (.h5) format. Default saves output in current dir.

    Args:
        filename (str): Name of filterbank file to read
        out_dir (str): Output directory path. Defaults to cwd
        new_filename (None or str): Name of output filename. If not set, will default
                                    to same as input, but with .h5 instead of .fil
        t_start (int): Start integration ID to be extracted from file
        t_stop (int): Stop integration ID to be extracted from file
    """

    wf = Waterfall(filename, load_data=False, t_start=t_start, t_stop=t_stop)
    if not new_filename:
        new_filename = out_dir+filename.replace('.fil', '.h5').split('/')[-1]

    if '.h5' not in new_filename:
        new_filename = new_filename+'.h5'

    write_to_h5(wf, new_filename)


def cmd_tool(args=None):
    """ Command line utility for converting Sigproc filterbank (.fil) to  HDF5 (.h5) format

    Usage:
        fil2h5 <FULL_PATH_TO_FIL_FILE> [options]

    Options:
      -h, --help            show this help message and exit
      -o OUT_DIR, --out_dir=OUT_DIR
                            Location for output files. Default: local dir.
      -n NEW_FILENAME, --new_filename=NEW_FILENAME
                            New name. Default: replaces extention to .h5
      -d, --delete_input    This option deletes the input file after conversion.
      -l MAX_LOAD           Maximum data limit to load. Default:1GB
    """

    parser = ArgumentParser(description="Command line utility for converting Sigproc filterbank (.fil) to  HDF5 (.h5) format  \n >>fil2h5 <FULL_PATH_TO_FIL_FILE> [options]")
    parser.add_argument("filepath_in", type=str, help="Path of input Filterbank file")
    parser.add_argument('-o', '--out_dir', dest='out_dir', type=str, default='./', help='Location for output files. Default: local dir. ')
    parser.add_argument('-n', '--new_filename', dest='new_filename', type=str, default='', help='New name. Default: replaces extention to .h5')
    parser.add_argument('-d', '--delete_input', dest='delete_input', action='store_true', default=False, help='This option deletes the input file after conversion.')
    parser.add_argument('-s', '--start_id', dest='t_start', type=int, default=None, help='start integration ID')
    parser.add_argument('-t', '--stop_id', dest='t_stop', type=int, default=None, help='stop integration ID')

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    make_h5_file(args.filepath_in,
                 out_dir = args.out_dir,
                 new_filename = args.new_filename,
                 t_start=args.t_start,
                 t_stop=args.t_stop)

    if args.delete_input:
        logger.info("'Deleting input file: {}".format(args.filename_in))
        os.remove(args.filename_in)


if __name__ == "__main__":
    cmd_tool()
