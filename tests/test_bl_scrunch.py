"""
test_bl_scrunch
"""

import os
import pytest

import blimpy as bl
from tests.data import voyager_h5


def test_scrunch():
    """
    Tests the conversion of fil files into h5 in both light and heavy modes.
    But apparently it does not test the accuracy of the conversion.
    """

    # Creating test file.
    bl.bl_scrunch.bl_scrunch(voyager_h5, new_filename='test.scrunched.h5', f_scrunch=8)

    # Deleting test file
    os.remove('test.scrunched.h5')

def test_nameless():
    """
    The script should trigger a system
    exit if the user fails to provide a file name.
    """

    bl.bl_scrunch.bl_scrunch(voyager_h5, f_scrunch=8)
    
    with pytest.raises(SystemExit):
        bl.bl_scrunch.cmd_tool()

def test_cmd():
    """ Pass in some example sets of arguments """
    # this is a malformed arg set with one giant string in the first entry
    args = [voyager_h5 + ' -n \'test.scrunched.h5\' -d -f 8 -l .1']
    with pytest.raises(SystemExit):
        bl.bl_scrunch.cmd_tool(args)

if __name__ == "__main__":
    test_scrunch()
