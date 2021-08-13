r""" h5diag """

import os
import sys
from argparse import ArgumentParser
import h5py
from astropy.coordinates import Angle
from blimpy.io.hdf_reader import examine_h5


def oops(msg):
    print("\n*** h5diag: {}\n".format(msg))
    sys.exit(86)


def read_header(h5):
    """ Read header and return a Python dictionary of key:value pairs
    """

    header = {} # Initialise as a nil dictionary.

    for key, val in h5["data"].attrs.items():
        #if six.PY3:
        #    key = bytes(key, "ascii")
        if isinstance(val, bytes):
            val = val.decode("ascii")
        if key == "src_raj":
            header[key] = Angle(val, unit="hr")
        elif key == "src_dej":
            header[key] = Angle(val, unit="deg")
        else:
            header[key] = val

    return header


def examine(filename):
    r""" Diagnose the given HDF5 file"""
    h5 = h5py.File(filename, mode="r")
    examine_h5(h5)
    print("h5diag: data shape:", h5["data"].shape)


def cmd_tool(args=None):
    """ Command line tool h5diag """
    parser = ArgumentParser(description="Command line utility for diagnosing HDF5 files.")
    parser.add_argument("filename", type=str, help="Path of file to read")
    if args is None:
        args = sys.argv[1:]
    parse_args = parser.parse_args(args)

    if not os.path.isfile(parse_args.filename):
        oops("Not a file: {}".format(parse_args.filename))
    if not h5py.is_hdf5(parse_args.filename):
        oops("Not an HDF5 file: {}".format(parse_args.filename))

    print("h5diag: Begin")
    examine(parse_args.filename)
    print("h5diag: No errors detected")


if __name__ == "__main__":
    cmd_tool()
