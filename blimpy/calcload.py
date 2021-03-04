r''' calcload.py - Calculate the Waterfall max_load value needed to load the data array for a given file.'''

from os.path import getsize
import numpy as np
from argparse import ArgumentParser
import blimpy as bl


def calc_max_load(arg_path):
    r''' Calculate the max_load parameter value for a subsequent Waterfall instantiation.
    
    Algorithm:
        * A = file size.
        * B = data array size (blimpy currently makes a copy of the data array)
        * Return ceil(A + B in GB)
    '''
    wf = bl.Waterfall(arg_path, load_data=False)
    data_size = float(wf.header['nchans'] * wf.n_ints_in_file * wf.header['nbits']) / 8.0
    ngbytes = (float(getsize(arg_path)) + data_size) / 1e9
    max_load = np.ceil(ngbytes)
    print('plot_event calc_max_load: max_load={} is required for {}'.format(max_load, arg_path))
    return max_load


def cmd_tool(args=None):
    p = ArgumentParser(description='Calculate the Waterfall max_load value needed to load the data array for a given file.')
    p.add_argument('filepath', type=str, help='Name of filepath to open (h5 or fil)')
    
    if args is None:
        args = p.parse_args()
    else:
        args = p.parse_args(args)

    gb = calc_max_load(args.filepath)
    if gb > 1.0:
        print('Use Waterfall instantiation with a max_load={}'.format(gb))
    else:
        print('Use Waterfall without a max_load= specification')


if __name__ == "__main__":
    cmd_tool()