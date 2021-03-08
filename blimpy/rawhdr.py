r'''
Read the specified raw file.
Print the entire header.
Highlight the number of observation channels.
'''

import sys
import blimpy
from argparse import ArgumentParser


def get_obsnchan(filepath):
    gr = blimpy.GuppiRaw(filepath)
    header, data_idx = gr.read_header()
    return header['OBSNCHAN']


def cmd_tool(args=None):
    p = ArgumentParser(description='Calculate the Waterfall max_load value needed to load the data array for a given file.')
    p.add_argument('filepath', type=str, help='Name of raw guppi file path to access')   
    if args is None:
        args = p.parse_args()
    else:
        args = p.parse_args(args)

    print('Getting the header from {}.\n'.format(args.filepath))
    with open(args.filepath, 'rb') as fh:
        while True:
            buffer = fh.read(80).decode("utf-8").strip()

            print('\t', buffer)
            if buffer[0:3] == "END":
                break
    print('\nThe number of observation channels:', get_obsnchan(args.filepath))

if __name__ == "__main__":
    cmd_tool()
