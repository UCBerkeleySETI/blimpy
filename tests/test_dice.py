from blimpy import dice
import pytest


def test_dice():
    with pytest.raises(SystemExit):
        dice.cmd_tool()


if __name__ == "__main__":
    test_dice()