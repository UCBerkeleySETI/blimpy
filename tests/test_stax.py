from os.path import dirname
import pytest
from blimpy.stax import cmd_tool
from tests.data import voyager_fil, voyager_h5


def test_stax():
    dir = dirname(voyager_h5)
    args = [voyager_fil, voyager_h5, "--plot_dir", dir]
    cmd_tool(args)
    args = [voyager_fil, voyager_h5, "--plot_dir", dir, "--f_start", "8419", "--f_stop", "8420"]
    cmd_tool(args)


if __name__ == "__main__":
    test_stax()
