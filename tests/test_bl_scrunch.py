r"""
test_bl_scrunch
"""

import os
import pytest

import blimpy as bl
from tests.data import voyager_h5


def test_scrunch():
    r"""
    Tests the conversion of fil files into h5 in both light and heavy modes.
    But apparently it does not test the accuracy of the conversion.
    """
    print("\n===== test_scrunch BEGIN")
    # Creating test file.
    bl.bl_scrunch.bl_scrunch(voyager_h5, new_filename='test.scrunched.h5', f_scrunch=8)

    # Deleting test file
    os.remove('test.scrunched.h5')
    print("\n===== test_scrunch END")


def test_nameless():
    r"""
    The script should trigger a system
    exit if the user fails to provide a file name.
    """

    print("\n===== test_nameless BEGIN")
    bl.bl_scrunch.bl_scrunch(voyager_h5, f_scrunch=8)

    with pytest.raises(SystemExit):
        bl.bl_scrunch.cmd_tool()
    print("\n===== test_nameless END")


def test_cmd():
    r""" Pass in some example sets of arguments """
    # this is a malformed arg set with one giant string in the first entry
    print("\n===== test_cmd BEGIN")
    args = [voyager_h5, "-n", "\'test.scrunched.h5\'", "-f", "8", "-l", "0.1"]
    bl.bl_scrunch.cmd_tool(args)
    args = [voyager_h5, "-f", "8", "-l", "0.1"]
    bl.bl_scrunch.cmd_tool(args)
    args = [voyager_h5, "--mickey_mouse", "-f", "8", "-l", "0.1"]
    with pytest.raises(SystemExit):
        bl.bl_scrunch.cmd_tool(args)
    print("\n===== test_cmd END")


if __name__ == "__main__":
    test_scrunch()
    test_nameless()
    test_cmd()
