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

    # Testing filename
    bl.fil2h5.make_h5_file(voyager_fil, new_filename='test')

    # Deleting test file
    os.remove('test.h5')

def test_cmd_tool():
    """
    This is the same test file, but now through the cmd tool.
    """
    #with pytest.raises(SystemExit):
    bl.fil2h5.cmd_tool([voyager_fil, '-n', 'cmd.h5'])

def test_no_args():
    """
    The cmd tool needs to exit, mandating a file name.
    """
    with pytest.raises(SystemExit):
        bl.fil2h5.cmd_tool()


if __name__ == "__main__":
    test_fil2h5_conversion()
