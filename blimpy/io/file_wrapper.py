#!/usr/bin/env python
""" This model handles file types.
"""

import os
import sys
import h5py
import six

# This import relates to a failing project, so we will
# probably want to remove the import at some point in the future.
import blimpy.io.sigproc

# import pdb;# pdb.set_trace()

from blimpy.io.fil_reader import FilReader
from blimpy.io.hdf_reader import H5Reader

def open_file(filename, f_start=None, f_stop=None,t_start=None, t_stop=None,load_data=True,max_load=1.):
    """Open a HDF5 or filterbank file

    Returns instance of a Reader to read data from file.

    ================== ==================================================
    Filename extension File type
    ================== ==================================================
    h5, hdf5           HDF5 format
    fil                fil format
    *other*            Will raise NotImplementedError
    ================== ==================================================

    """
    if not os.path.isfile(filename):
        type(filename)
        print(filename)
        raise IOError("No such file or directory: " + filename)

    filename = os.path.expandvars(os.path.expanduser(filename))
    # Get file extension to determine type
    ext = filename.split(".")[-1].strip().lower()

    if six.PY3:
        ext = bytes(ext, 'ascii')

    if h5py.is_hdf5(filename):
        # Open HDF5 file
        return H5Reader(filename, f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop,
                        load_data=load_data, max_load=max_load)
    elif blimpy.io.sigproc.is_filterbank(filename):
        # Open FIL file
        return FilReader(filename, f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop, load_data=load_data, max_load=max_load)
    else:
        raise NotImplementedError('Cannot open this type of file with Waterfall')
