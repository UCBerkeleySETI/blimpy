import os
from blimpy import dice
import pytest

from tests.data import voyager_h5
from tests.data import voyager_fil
TEST_DATA_DIR = os.path.dirname(voyager_h5)
HERE = os.getcwd()

# This section makes sure that the tool exits correctly
# when certain issues are encountered.

def test_no_args():
    """
    Make sure the tool closes if no arguments are passed in
    """
    with pytest.raises(SystemExit):
        dice.cmd_tool()

def test_no_input_file():
    with pytest.raises(SystemExit):
        args = ['-b', '8419.24', '-e', '8419.35', '-x', 'h5', '-o', 'test_dice.h5']
        dice.cmd_tool(args)

def test_missing_format_type():
    with pytest.raises(SystemExit):
        args = ['-f', voyager_h5, '-b', '8419.24', '-e', '8419.35', '-x']
        dice.cmd_tool(args)

# This section makes sure that the tool can handle various
# file formats without exception.

def test_h5():
    os.chdir(TEST_DATA_DIR)
    args = ['-f', voyager_h5, '-b', '8419.24', '-e', '8419.35', '-x', 'h5', '-o', 'test_dice.h5']
    dice.cmd_tool(args)
    os.chdir(HERE)

def test_h5_no_out_file():
    os.chdir(TEST_DATA_DIR)
    args = ['-f', voyager_h5, '-b', '8419.24', '-e', '8419.35', '-x', 'h5']
    dice.cmd_tool(args)
    os.chdir(HERE)

def test_fil():
    os.chdir(TEST_DATA_DIR)
    args = ['-f', voyager_fil, '-b', '8419.24', '-e', '8419.35', '-x', 'fil', '-o', 'test_dice.fil']
    dice.cmd_tool(args)
    os.chdir(HERE)

def test_fil_no_out_file():
    os.chdir(TEST_DATA_DIR)
    args = ['-f', voyager_fil, '-b', '8419.24', '-e', '8419.35', '-x', 'fil']
    dice.cmd_tool(args)
    os.chdir(HERE)

