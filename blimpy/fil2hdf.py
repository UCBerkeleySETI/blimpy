#!/usr/bin/env python

from .filterbank import Filterbank
import h5py
try:
    HAS_BITSHUFFLE = True
    import bitshuffle.h5
except ImportError:
    HAS_BITSHUFFLE = False
import time
import os
import glob

import numpy as np


MAX_SIZE = 4e9

def cmd_tool(args=None):
    """ Command line utility for creating HDF5 blimpy files. """
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Command line utility for creating HDF5 Filterbank files.")
    parser.add_argument('dirname', type=str, help='Name of directory to read')
    args = parser.parse_args()
    
    if not HAS_BITSHUFFLE:
        print("Error: the bitshuffle library is required to run this script.")
        exit()

    filelist = glob.glob(os.path.join(args.dirname, '*.fil'))

    for filename in filelist:
        if not os.path.exists(filename + '.h5'):
            t0 = time.time()
            print("\nReading %s header..." % filename)
            fb = Filterbank(filename, load_data=False)

            data_shape = (fb.n_ints_in_file, fb.header['nifs'], fb.header['nchans'])
            data_dtype = fb.data.dtype
            print(data_dtype)


            print("Creating new dataset, %s" % str(data_shape))
            block_size = 0
            h5 = h5py.File(filename + '.h5', 'w')

            h5.attrs['CLASS'] = 'FILTERBANK'

            dset = h5.create_dataset('data',
                              shape=data_shape,
                              compression=bitshuffle.h5.H5FILTER,
                              compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
                              dtype=data_dtype)

            dset_mask = h5.create_dataset('mask',
                                     shape=data_shape,
                                     compression=bitshuffle.h5.H5FILTER,
                                     compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
                                     dtype='uint8')

            dset.dims[0].label = "frequency"
            dset.dims[1].label = "feed_id"
            dset.dims[2].label = "time"

            dset_mask.dims[0].label = "frequency"
            dset_mask.dims[1].label = "feed_id"
            dset_mask.dims[2].label = "time"

            # Copy over header information as attributes
            for key, value in fb.header.items():
                dset.attrs[key] = value

            filesize = os.path.getsize(filename)

            if filesize >= MAX_SIZE:
                n_int_per_read = int(filesize / MAX_SIZE / 2)
                print("Filling in with data over %i reads..." % n_int_per_read)
                for ii in range(0, n_int_per_read):
                    print("Reading %i of %i" % (ii + 1, n_int_per_read))
                    #print  ii*n_int_per_read, (ii+1)*n_int_per_read
                    fb = Filterbank(filename, t_start=ii*n_int_per_read, t_stop=(ii+1) * n_int_per_read)
                    dset[ii*n_int_per_read:(ii+1)*n_int_per_read] = fb.data[:]
            else:
                fb = Filterbank(filename)
                print(dset.shape, " -> ", fb.data.shape)
                dset[:] = fb.data[:]

            h5.close()

            t1 = time.time()
            print("Conversion time: %2.2fs" % (t1- t0))

if __name__ == "__main__":
    cmd_tool()