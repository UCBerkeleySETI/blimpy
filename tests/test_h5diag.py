from os.path import dirname
import numpy as np
import hdf5plugin
import h5py
from blimpy.h5diag import cmd_tool
from tests.data import voyager_h5, voyager_fil
import pytest


header = [
    ["fruit", "apple"],
    ["color", "red"],
    ["plant", "tree"]
    ]
DIR = dirname(voyager_fil)
TEST_H5 = DIR + "/test.h5"
TIME_INSTANCES = 8
FREQ_INSTANCES = 16
DATA_BYTESIZE = TIME_INSTANCES * FREQ_INSTANCES * 4


def my_writer(my_class):
    data_out = np.ndarray(shape=(TIME_INSTANCES, 1, FREQ_INSTANCES), dtype=float)
    for ii in range(TIME_INSTANCES):
        for jj in range(FREQ_INSTANCES):
            data_out[ii, 0, jj] = 42.0
    print("data_out shape:", data_out.shape)
    
    with h5py.File(TEST_H5, "w") as h5:
        h5.attrs["CLASS"]   = my_class
        h5.attrs["VERSION"] = "1.0"
    
        bs_compression = hdf5plugin.Bitshuffle(nelems=0, lz4=True)["compression"]
        bs_compression_opts = hdf5plugin.Bitshuffle(nelems=0, lz4=True)["compression_opts"]
    
        dset = h5.create_dataset("data",
                                 data=data_out,
                                 compression=bs_compression,
                                 compression_opts=bs_compression_opts)
    
        dset_mask = h5.create_dataset("mask",
                                      shape=data_out.shape,
                                      compression=bs_compression,
                                      compression_opts=bs_compression_opts,
                                      dtype="uint8")
    
        dset.dims[2].label = b"frequency"
        dset.dims[1].label = b"feed_id"
        dset.dims[0].label = b"time"
    
        dset_mask.dims[2].label = b"frequency"
        dset_mask.dims[1].label = b"feed_id"
        dset_mask.dims[0].label = b"time"
    
        # Copy over header information as attributes
        for key, value in header:
            dset.attrs[key] = value


def execute_command(args):
    print("\ntest_h5diag: args:", args)
    cmd_tool(args)


def test_h5diag():

    args = [voyager_h5]
    execute_command(args)

    with pytest.raises(SystemExit):
        args = [voyager_fil]
        execute_command(args)
    
    my_writer("FRUITY")
    with pytest.raises(SystemExit):
        args = [TEST_H5]
        execute_command(args)
    
    with h5py.File(TEST_H5, "w") as h5:
        h5.attrs["VERSION"] = "42.0"
    with pytest.raises(SystemExit):
        args = [TEST_H5]
        execute_command(args)
    
    with h5py.File(TEST_H5, "w") as h5:
        h5.attrs["CLASS"] = "FILTERBANK"
    with pytest.raises(SystemExit):
        args = [TEST_H5]
        execute_command(args)

    with h5py.File(TEST_H5, "w") as h5:
        h5.attrs["CLASS"] = "FILTERBANK"
        h5.attrs["VERSION"] = "42.0"
    with pytest.raises(SystemExit):
        args = [TEST_H5]
        execute_command(args)

    my_writer("FILTERBANK")
    args = [TEST_H5]
    execute_command(args)


if __name__ == "__main__":
    test_h5diag()
