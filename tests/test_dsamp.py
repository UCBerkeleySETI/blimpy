"""
# test_dsamp
"""

import pytest
import blimpy as bl
from tests.data import voyager_fil, voyager_h5, test_h5


GROUP_SIZE = 3


def test_dsamp_fil_to_h5():
    """ fil to h5 test.
    """
    bl.dsamp.make_h5_file(voyager_fil, test_h5, GROUP_SIZE)


def test_dsamp_h5_to_h5():
    """ h5 to h5 test.
    """
    bl.dsamp.make_h5_file(voyager_h5, test_h5, GROUP_SIZE)


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
