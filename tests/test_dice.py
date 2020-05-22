from blimpy import dice
import pytest

from tests.data import voyager_h5
from tests.data import voyager_fil

# This section makes sure that the tool exits correctly
# when certain issues are encountered.

def test_no_args():
    """
    Make sure the tool closes if no arguments are passed in
    """
    with pytest.raises(SystemExit):
        dice.cmd_tool()

def test_no_in_file():
    """
    This test does not actually increase coverage,
    which causes me to doubt the effectiveness
    of this approach to testing for errors...
    """
    with pytest.raises(SystemExit):
        args = ['-b 8419.24', '-e 8419.35', '-x', 'h5', '-otest_dice.h5']
        dice.cmd_tool()
    with pytest.raises(SystemExit):
        args = ['-f%s' % voyager_h5, '-b 8419.24', '-e 8419.35', '-x']
        dice.cmd_tool()

# This section makes sure that the tool can handle various
# file formats without exception.

def test_h5():
    args = ['-f%s' % voyager_h5, '-b 8419.24', '-e 8419.35', '-x', 'h5', '-otest_dice.h5']
    dice.cmd_tool(args)

def test_h5_no_out_file():
    args = ['-f%s' % voyager_h5, '-b 8419.24', '-e 8419.35', '-x', 'h5']
    dice.cmd_tool(args)

def test_fil():
    args = ['-f%s' % voyager_fil, '-b 8419.24', '-e 8419.35', '-x', 'fil', '-otest_dice.fil']
    dice.cmd_tool(args)

def test_fil_no_out_file():
    args = ['-f%s' % voyager_fil, '-b 8419.24', '-e 8419.35', '-x', 'fil']
    dice.cmd_tool(args)

if __name__ == "__main__":
    test_dice()
