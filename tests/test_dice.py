from blimpy import dice
import pytest

from tests.data import voyager_h5

def test_dice():
    with pytest.raises(SystemExit):
        dice.cmd_tool()

    args = ['-f%s' % voyager_h5, '-b 8419.24', '-e 8419.35', '-x', 'h5', '-otest_dice.h5']
    dice.cmd_tool(args)

if __name__ == "__main__":
    test_dice()