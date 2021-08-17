import pytest
import blimpy as bl
from tests.data import voyager_fil

def test_reset_works():
    """
    The test is trivial, but coverage must be
    increased at all costs.
    """
    a, b = bl.match_fils.reset_outs()
    assert a is None
    assert b is None

def test_batcher():
    """
    This test may be obsolete once/if test_cmd
        can be patched up.
    """
    bl.match_fils.make_batch_script()

def test_find_header_size():
    assert bl.match_fils.find_header_size(voyager_fil) > 0

def test_cmdline():
    args = [voyager_fil, voyager_fil]
    with pytest.raises(OSError):
        bl.match_fils.cmd_tool(args)

