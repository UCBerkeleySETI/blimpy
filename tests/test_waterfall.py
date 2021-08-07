import os
import numpy as np
from tests.data import voyager_h5, voyager_fil
import blimpy as bl

import pytest

OUTDIR = os.path.dirname(voyager_h5) + "/"

def test_info():
    a = bl.Waterfall(voyager_h5)
    print(a)
    a.info()
    a.blank_dc(n_coarse_chan=1)
    a.calibrate_band_pass_N1()
    # Below: orphaned functions (not used anywhere else).
    # That is a little strange...
    a.grab_data()
    a.read_data()
    # It is VERY strange for INTERNAL functions to be orphaned!
    a._get_chunk_dimensions()
    # plenty of missed if suites
    a._get_blob_dimensions((300, 300, 300, 300))
    a._update_header()

def test_get_freqs():
    wf = bl.Waterfall(voyager_h5)
    freqs = wf.container.populate_freqs()
    sum1 = np.sum(freqs)
    freqs = wf.get_freqs()
    sum2 = np.sum(freqs)
    assert sum1 == sum2
    wf = bl.Waterfall(voyager_fil)
    freqs = wf.container.populate_freqs()
    sum1 = np.sum(freqs)
    freqs = wf.get_freqs()
    sum2 = np.sum(freqs)
    assert sum1 == sum2

def test_cmdline():
    from blimpy.waterfall import cmd_tool

    args = [voyager_h5, '-S', '-s', OUTDIR + 'test.png']
    cmd_tool(args)

    args = [voyager_h5, '-i']
    cmd_tool(args)

    if os.path.exists(OUTDIR + 'test.h5'):
        os.remove(OUTDIR + 'test.h5')
        args = [voyager_h5, '-H', '-o', OUTDIR + 'test.h5']
        cmd_tool(args)
        assert os.path.exists(OUTDIR + 'test.h5')
        os.remove(OUTDIR + 'test.h5')

    if os.path.exists(OUTDIR + 'test.fil'):
        os.remove(OUTDIR + 'test.fil')
        args = [voyager_h5, '-F', '-o', OUTDIR + 'test.fil']
        cmd_tool(args)
        assert os.path.exists(OUTDIR + 'test.fil')
        os.remove(OUTDIR + 'test.fil')

def test_cmd_arguments():
    from blimpy.waterfall import cmd_tool
    
    args = [voyager_h5, '-H', '-F', '-o', OUTDIR + 'test.fil']
    with pytest.raises(ValueError):
        cmd_tool(args)

if __name__ == "__main__":
    test_info()
    test_cmdline()
    test_cmd_arguments()
