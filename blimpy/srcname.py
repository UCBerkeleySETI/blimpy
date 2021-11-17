""" srcname 

    Change the HDF5 Filterbank header source_name field value.
"""

import os
import sys
from argparse import ArgumentParser
import h5py
from astropy.coordinates import Angle


def oops(msg):
    print("\n*** srcname: {}\n".format(msg))
    sys.exit(86)


def read_header(h5):
    """
    Read the Filterbank header and return a Python dictionary of key:value pairs

    Parameters
    ----------
    h5 : HDF5 file handle
        This represents an open HDF5 file.

    Returns
    -------
    header : dict
        The Filterbank header.

    """
    header = {} # Initialise as a nil dictionary.

    for key, val in h5["data"].attrs.items():
        if isinstance(val, bytes):
            val = val.decode("ascii")
        if key == "src_raj":
            header[key] = Angle(val, unit="hr")
        elif key == "src_dej":
            header[key] = Angle(val, unit="deg")
        else:
            header[key] = val

    return header


def examine(filepath):
    """
    Diagnose the specified HDF5 file.

    Parameters
    ----------
    filepath : path
        An O/S path to an HDF5 file.

    Returns
    -------
    header : dict
        The Filterbank header.

    """
    h5 = h5py.File(filepath, mode="r")
    if "CLASS" in h5.attrs:
        classstr = h5.attrs["CLASS"]
    else:
        oops("CLASS attribute missing")
    if classstr != "FILTERBANK":
        oops("Expected CLASS attribute to be 'FILTERBANK' but saw '{}'".format(classstr))
    if "VERSION" in h5.attrs:
        versionstr = h5.attrs["VERSION"]
    else:
        oops("VERSION attribute missing")
    header = read_header(h5)
    print("Header:", header)
    if not "data" in h5:
        oops("data attribute missing")
    if h5["data"].ndim != 3:
        oops("Expected data.ndim to be 3 but saw '{}'".format(h5["data"].ndim))
    print("srcname: data shape:", h5["data"].shape)
    h5.close()

    return header


def cmd_tool(args=None):
    """ Command line tool srcname """
    parser = ArgumentParser(description="Patch the header source field in an HDF5 file.")
    parser.add_argument("filepath", default=None, type=str, help="Path of file to read")
    parser.add_argument("new_source_name", default=None, type=str, 
                        help="New header source name field value")
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    if not os.path.isfile(args.filepath):
        oops("Not a file: {}".format(args.filepath))
    if not h5py.is_hdf5(args.filepath):
        oops("Not an HDF5 file: {}".format(args.filepath))

    header = examine(args.filepath)
    print("srcname: No errors detected in {}".format(args.filepath))
    print("\nThe current source_name field is [{}]".format(header["source_name"]))
    input("Are you sure you want to replace it with [{}]?  Press Enter to continue.  Ctrl-C to cancel: "
          .format(args.new_source_name))
    h5 = h5py.File(args.filepath, mode="r+")
    h5["data"].attrs["source_name"] = args.new_source_name
    h5.close()
    print("All done, best wishes!")


if __name__ == "__main__":
    cmd_tool()
