"""
Test the calcload.py utility and function calc_max_load()
"""
from tests.data import voyager_h5, voyager_fil
from blimpy.calcload import cmd_tool, calc_max_load

def test_calcload():
    r""" Test the calcload command line tool """
    args = [voyager_h5]
    cmd_tool(args)
    args = ['-v', voyager_fil]
    cmd_tool(args)

def test_calc_max_load():
    gb1 = calc_max_load(voyager_h5)
    gb2 = calc_max_load(voyager_fil)
    assert(gb1 == gb2 == 1.0)
   
if __name__ == "__main__":
    test_calcload()
    test_calc_max_load()
