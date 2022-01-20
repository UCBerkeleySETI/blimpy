"""
# test_h52fil
"""

import pytest

import os
import blimpy as bl
from tests.data import voyager_h5


VOYA_DIR = os.path.dirname(voyager_h5) + "/"


def name_case(in_string, out_string):
    infile = in_string
    outfile = out_string
    os.system("cp " + voyager_h5 + " " + infile)
    bl.h52fil.make_fil_file(infile)
    if not os.path.exists(outfile):
        print("\n*** name_case: file {} does not exist.  Input file {}\n".format(outfile, infile))
        assert False
    os.remove(infile)
    os.remove(outfile)

def test_h52fil_conversion():
    """ Tests the conversion of fil files into h5 in both light and heavy modes.
    """

    # Creating test file.
    bl.h52fil.make_fil_file(voyager_h5, new_filename='test.fil')

    # Creating a "large" test file.
    bl.h52fil.make_fil_file(voyager_h5, new_filename='test_large.fil', max_load=0.001)

    # Testing filename
    bl.h52fil.make_fil_file(voyager_h5, new_filename='test')

    # Deleting test file
    os.remove('test.fil')
    os.remove('test_large.fil')

def test_help():
    """
    The user of these tests should verify that the help statement
    was printed as expected; this test merely verifies that the
    system exited.
    """
    with pytest.raises(SystemExit):
        bl.h52fil.cmd_tool(['-h'])

def test_cmd_tool():
    """
    This is the same Voyager test file, but now through the cmd tool.
    """
    bl.h52fil.cmd_tool([voyager_h5, '-n', VOYA_DIR + 'cmd.fil', '-l', '.001'])

def test_no_args():
    """
    The cmd tool needs to exit, mandating a file name.
    """
    with pytest.raises(SystemExit):
        bl.h52fil.cmd_tool(['-n', VOYA_DIR + 'cmd.fil', '-l', '.001'])

def test_fil2h5_input_names():
    """ Make sure that the output name does not get mangled.
    """
    name_case("abcd.filter.def.h5", "abcd.filter.def.fil")
    name_case("abcd.efgh", "abcd.efgh.fil")
    name_case("abcd", "abcd.fil")

if __name__ == "__main__":
    test_h52fil_conversion()
    test_help()
    test_cmd_tool()
    test_no_args()
    test_fil2h5_input_names()
    
