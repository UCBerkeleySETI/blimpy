"""
# test_fil2h5
"""

import pytest

import os
import blimpy as bl
from tests.data import voyager_fil


def test_fil2h5_conversion():
    """ Tests the conversion of fil files into h5 in both light and heavy modes.
    """

    # Creating test file.
    bl.fil2h5.make_h5_file(voyager_fil, new_filename='test.h5')

    # Creating a "large" test file.
    bl.fil2h5.make_h5_file(voyager_fil, new_filename='test_large.h5', max_load=0.001)

    # Testing filename
    bl.fil2h5.make_h5_file(voyager_fil, new_filename='test')

    # Deleting test file
    os.remove('test.h5')
    os.remove('test_large.h5')

def test_help():
    """
    The user of these tests should verify that the help statement
    was printed as expected; this test merely verifies that the
    system exited.
    """
    with pytest.raises(SystemExit):
        bl.fil2h5.cmd_tool(['-h'])

def test_cmd_tool():
    """
    This is the same large test file, but now through the cmd tool.
    """
    #with pytest.raises(SystemExit):
    bl.fil2h5.cmd_tool([voyager_fil, '-n', 'cmd.h5', '-l', '0.001'])

def test_no_args():
    """
    The cmd tool needs to exit, mandating a file name.
    """
    with pytest.raises(SystemExit):
        bl.fil2h5.cmd_tool()


if __name__ == "__main__":
    test_fil2h5_conversion()
