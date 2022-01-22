"""
# test_fil2h5
"""

import pytest

import os
import blimpy as bl
from tests.data import voyager_fil


VOYA_DIR = os.path.dirname(voyager_fil) + "/"

def name_case(in_string, out_string):
    infile = in_string
    outfile = out_string
    os.system("cp " + voyager_fil + " " + infile)
    bl.fil2h5.make_h5_file(infile)
    if not os.path.exists(outfile):
        print("\n*** name_case: file {} does not exist.  Input file {}\n".format(outfile, infile))
        assert False
    os.remove(infile)
    os.remove(outfile)

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
    args = [voyager_fil, '-n', VOYA_DIR + 'cmd.h5']
    bl.fil2h5.cmd_tool(args=args)

def test_fil2h5_input_names():
    """ Make sure that the output name does not get mangled.
    """
    name_case("abcd.filter.def.fil", "abcd.filter.def.h5")
    name_case("abcd.efgh", "abcd.efgh.h5")
    name_case("abcd", "abcd.h5")

def test_no_args():
    """
    The cmd tool needs to exit, mandating a file name.
    """
    with pytest.raises(SystemExit):
        bl.fil2h5.cmd_tool("")


if __name__ == "__main__":
    test_fil2h5_conversion()
    test_cmd_tool()
    test_fil2h5_input_names()
    test_no_args()
