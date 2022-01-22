''' calcload.py - Calculate the Waterfall max_load value needed to load the data array for a given file.'''

import sys
from argparse import ArgumentParser
import numpy as np
import blimpy as bl


def calc_max_load(arg_path, verbose=False):
    r''' Calculate the max_load parameter value for a subsequent Waterfall instantiation.

    Algorithm:
        * A = minimum Waterfall object size.
        * B = data array size within one polarisation.
        * Return ceil(A + B in GB)
    '''
    wf = bl.Waterfall(arg_path, load_data=False)
    min_size = float(sys.getsizeof(wf.header)) + float(sys.getsizeof(wf))
    data_size = float(wf.header['nchans'] * wf.n_ints_in_file * wf.header['nbits']) / 8.0
    ngbytes = (min_size + data_size) / 1e9
    max_load = np.ceil(ngbytes)
    if verbose:
        print('calc_max_load: Waterfall object size excluding data = {}, data array size = {}, total GBs = {:.1f}'
              .format(min_size, data_size, ngbytes))
    return max_load


def cmd_tool(args=None):
    r'''Command line entrypoint for "calcload"'''
    p = ArgumentParser(description='Calculate the Waterfall max_load value needed to load the data array for a given file.')
    p.add_argument('filepath', type=str, help='Name of filepath to open (h5 or fil)')
    p.add_argument('-v', action='store_true', default=False, dest='verbose',
                       help='verbose output if True.')

    if args is None:
        args = p.parse_args()
    else:
        args = p.parse_args(args)

    gb = calc_max_load(args.filepath, args.verbose)
    if gb > 1.0:
        print('Use Waterfall instantiation with a max_load={}'.format(gb))
    else:
        print('Use Waterfall without a max_load= specification')


if __name__ == "__main__":
    cmd_tool()
