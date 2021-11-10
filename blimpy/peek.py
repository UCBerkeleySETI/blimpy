"""
Simple script for quickly peeking at the data matrix.
"""

import sys
from argparse import ArgumentParser

from blimpy import Waterfall


def oops(msg):
    print("\n*** OOPS, {} !!!".format(msg))
    sys.exit(86)


def cmd_tool(args=None):
    """ Command line utility for peeking at the data matrix.
    """

    parser = ArgumentParser(description="Command line utility for peeking at a .fil or .h5 file")
    parser.add_argument("filepath_in", type=str, help="Path of input .fil or .h5 file")
    parser.add_argument("-i", "--start_tint", dest="start_tint", type=int, default=0,
                        help="Starting integration index relative to 0.  Default: 0")
    parser.add_argument("-Z", "--if", dest="the_IF", type=int, default=0,
                        help="Starting IF index relative to 0.  Default: 0")
    parser.add_argument("-c", "--start_fchan", dest="start_fchan", type=int, default=0,
                        help="Starting fine channel index relative to 0.  Default: 0")


    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    wf = Waterfall(args.filepath_in)

    if args.start_tint < 0 or args.start_tint > (wf.n_ints_in_file - 2):
        oops("--start_tint is not a valid time integration index")
    if args.the_IF < 0:
        oops("--if is not a valid IF index")
    if args.start_fchan < 0 or args.start_fchan > (wf.header["nchans"] - 3):
        oops("--start_fchan is not a valid fine channel index")

    print("Fine channel frequency columns     go this way ------->")
    print("Integration_{}:                     {}   {}   {}"
          .format(args.start_tint,
                  wf.data[args.start_tint, args.the_IF, args.start_fchan],
                  wf.data[args.start_tint, args.the_IF, args.start_fchan + 1],
                  wf.data[args.start_tint, args.the_IF, args.start_fchan + 2]))
    print("Integration_{}:                     {}   {}   {}"
          .format(args.start_tint + 1,
                  wf.data[args.start_tint + 1, args.the_IF, args.start_fchan],
                  wf.data[args.start_tint + 1, args.the_IF, args.start_fchan + 1],
                  wf.data[args.start_tint + 1, args.the_IF, args.start_fchan + 2]))


if __name__ == "__main__":
    cmd_tool()
