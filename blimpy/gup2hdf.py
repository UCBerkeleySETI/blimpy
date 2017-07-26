#!/usr/bin/env python

from .guppi import GuppiRaw
import h5py
try:
    import bitshuffle.h5
    HAS_BITSHUFFLE = True
except ImportError:
    HAS_BITSHUFFLE = False
    
import time
import os
import glob
import numpy as np

def cmd_tool(args=None):
    """ Command line tool for converting guppi raw into HDF5 versions of guppi raw """
    from argparse import ArgumentParser

    if not HAS_BITSHUFFLE:
        print("Error: the bitshuffle library is required to run this script.")
        exit()

    parser = ArgumentParser(description="Command line utility for creating HDF5 Raw files.")
    parser.add_argument('filename', type=str, help='Name of filename to read')
    args = parser.parse_args()

    fileroot = args.filename.split('.0000.raw')[0]

    filelist = glob.glob(fileroot + '*.raw')
    filelist = sorted(filelist)


    # Read first file
    r = GuppiRaw(filelist[0])
    header, data = r.read_next_data_block()
    dshape = data.shape #r.read_next_data_block_shape()
    print(dshape)

    n_blocks_total = 0
    for filename in filelist:
        print(filename)
        r = GuppiRaw(filename)
        n_blocks_total += r.n_blocks
    print(n_blocks_total)

    full_dshape = np.concatenate(((n_blocks_total,), dshape))


    # Create h5py file
    h5 = h5py.File(fileroot + '.h5', 'w')
    h5.attrs['CLASS'] = 'GUPPIRAW'
    block_size = 0      # This is chunk block size
    dset = h5.create_dataset('data',
                  shape=full_dshape,
                  #compression=bitshuffle.h5.H5FILTER,
                  #compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
                  dtype=data.dtype) 

    h5_idx = 0
    for filename in filelist:
        print("\nReading %s header..." % filename)
        r = GuppiRaw(filename)
        h5 = h5py.File(filename + '.h5', 'w')
    
        header, data = r.read_next_data_block()
        
        for ii in range(0, r.n_blocks):
            t0 = time.time()
            print("Reading block %i of %i" % (h5_idx+1, full_dshape[0]))
            header, data = r.read_next_data_block()
            t1 = time.time()
        
            t2 = time.time()
            print("Writing block %i of %i" % (h5_idx+1, full_dshape[0]))
            dset[h5_idx, :] = data
            t3 = time.time()
            print("Read: %2.2fs, Write %2.2fs" % ((t1-t0), (t3-t2)))
        
            h5_idx += 1

            # Copy over header information as attributes
            for key, value in header.items():
                dset.attrs[key] = value

        h5.close()

        t1 = time.time()
        print("Conversion time: %2.2fs" % (t1- t0))

if __name__ == "__main__":
    cmd_tool()