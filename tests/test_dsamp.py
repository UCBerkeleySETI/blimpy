"""
# test_dsamp
"""

import h5py
import numpy as np
import pytest
import blimpy as bl
from tests.data import voyager_fil, voyager_h5, test_h5, test_fil


GROUP_SIZE = 3


def check_h5_is_downsampled(in_path):
    """ The output .h5 data must be the sums of each group of input time samples.
    """
    in_data = bl.Waterfall(in_path, max_load=None).data
    n_out = in_data.shape[0] // GROUP_SIZE
    expected = in_data[:n_out * GROUP_SIZE].reshape(
        (n_out, GROUP_SIZE) + in_data.shape[1:]).sum(axis=1)
    with h5py.File(test_h5, 'r') as h5:
        assert np.allclose(h5['data'][:], expected)


def test_dsamp_fil_to_h5():
    """ fil to h5 test.
    """
    bl.dsamp.make_output_file(voyager_fil, test_h5, GROUP_SIZE, True)
    check_h5_is_downsampled(voyager_fil)


def test_dsamp_h5_to_h5():
    """ h5 to h5 test.
    """
    bl.dsamp.make_output_file(voyager_h5, test_h5, GROUP_SIZE, True)
    check_h5_is_downsampled(voyager_h5)


def test_dsamp_h5_to_fil():
    """ h5 to fil test.
    """
    bl.dsamp.make_output_file(voyager_h5, test_fil, GROUP_SIZE, False)


def test_dsamp_fil_to_fil():
    """ fil to fil test.
    """
    bl.dsamp.make_output_file(voyager_fil, test_fil, GROUP_SIZE, False)


def test_cmd_tool():
    """
    Exercise cmd_tool.
    """
    args = [voyager_fil, test_h5, "-s", str(GROUP_SIZE)]
    bl.dsamp.cmd_tool(args=args)


def test_no_args():
    """
    The cmd tool needs to exit, mandating a file name.
    """
    with pytest.raises(SystemExit):
        bl.dsamp.cmd_tool("")
    with pytest.raises(SystemExit):
        bl.dsamp.cmd_tool("-h")
