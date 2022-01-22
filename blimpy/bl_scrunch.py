r"""
    From an input HDF file or Filterbank file, perform channel averaging,
    producing a new HDF5 file of Filterbank file.
"""

import os
from argparse import ArgumentParser
from blimpy.waterfall import Waterfall
from .utils import change_the_ext


def bl_scrunch(in_path, out_dir='./', new_filename='', max_load=None, f_scrunch=None):
    r""" Frequency scrunch (lower resolution by averaging)

    Args:
        in_path : str
            Path of input file to open.
        out_dir : str
            Output directory.
        new_filename : str
            Output file name.
        max_load : int
            Waterfall object instantiation max_load parameter value.
        f_scrunch : int
            Number of frequency channels to average together at one time.
    """

    print("bl_scrunch: Input path: {}".format(in_path))
    in_ext = os.path.splitext(in_path)[1]
    if in_ext not in ('.fil', '.h5'):
        raise ValueError('Oops, input file extension must be .fil or .h5; saw: {} !'.format(in_ext))
    print("bl_scrunch: Averaging {} frequency channels at a time.".format(f_scrunch))

    wf = Waterfall(in_path, max_load=max_load)
    if new_filename == '':
        if in_ext == '.h5':
            out_path = out_dir + '/' + change_the_ext(in_path, 'h5', 'scrunched.h5').split('/')[-1]
        else: # .fil
            out_path = out_dir + '/' + change_the_ext(in_path, 'fil', 'scrunched.h5').split('/')[-1]
    else:
        out_path = out_dir + new_filename

    print("bl_scrunch: Output path: {}".format(out_path))
    wf.write_to_hdf5(out_path, f_scrunch=f_scrunch)
    print("bl_scrunch: End")


def cmd_tool(args=None):
    r"""  Command line utility for scrunching an input HDF5 file or Filterbank file.
    """

    p = ArgumentParser(description='Command line utility for converting HDF5 (.h5) to Sigproc filterbank (.fil) format \n >>h52fil <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_argument('filepath', type=str, help='Name of file path to open (.h5 or .fil).')
    p.add_argument('-f', '--fscrunch', dest='f_scrunch', type=int,
                 help='Number of frequency channels to average (scrunch) together.')
    p.add_argument('-o', '--out_dir', dest='out_dir', type=str, default='./',
                 help='Location for output files. Default: current directory.')
    p.add_argument('-n', '--new_filename', dest='new_filename', type=str, default='',
                 help='New filename. Default: replaces the file extension with .scrunched.fil or .scrunched .h5.')
    p.add_argument('-l', '--max_load', action='store', default=None, dest='max_load', type=float,
                 help='Maximum data limit to load. Default: 1.0 GB.')
    p.add_argument('-d', '--delete_input', dest='delete_input', action='store_true', default=False,
                 help='This option deletes the input file after conversion.')

    if args is None:
        args = p.parse_args()
    else:
        args = p.parse_args(args)

    bl_scrunch(args.filepath, out_dir=args.out_dir, new_filename=args.new_filename,
               max_load=args.max_load, f_scrunch=args.f_scrunch)

    if args.delete_input:
        print("'Deleting input file: %s"%(args.filepath))
        os.remove(args.filepath)


if __name__ == "__main__":

    cmd_tool()
