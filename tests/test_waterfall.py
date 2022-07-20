import os
import numpy as np
from tests.data import voyager_h5, voyager_fil
import blimpy as bl
from blimpy.waterfall import cmd_tool

import pytest

OUTDIR = os.path.dirname(voyager_h5) + "/"

def test_info():
    print("\n===== test_info")
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
    del a

def test_get_freqs():
    print("\n===== test_get_freqs")
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
    len_f = len(freqs)
    first = wf.header["fch1"]
    last_1 = freqs[-1]
    last_2 = first + (len_f - 1 ) *  wf.header["foff"]
    assert np.isclose(last_1, last_2, rtol=0.0001)

def test_cmdline():
    print("\n===== test_cmdline")

    args = [voyager_h5, '-S', '-p', 'w', '-s', OUTDIR + 'test.png']
    cmd_tool(args)

    args = [voyager_h5, '-S', '-p', 's', '-s', OUTDIR + 'test.png']
    cmd_tool(args)

    args = [voyager_h5, '-S', '-p', 'mm', '-s', OUTDIR + 'test.png']
    cmd_tool(args)

    args = [voyager_h5, '-S', '-p', 'k', '-s', OUTDIR + 'test.png']
    cmd_tool(args)

    args = [voyager_h5, '-S', '-p', 't', '-s', OUTDIR + 'test.png']
    cmd_tool(args)

    args = [voyager_h5, '-S', '-p', 'a', '-s', OUTDIR + 'test.png']
    cmd_tool(args)

    args = [voyager_h5, '-S', '-p', 'ank', '-s', OUTDIR + 'test.png']
    cmd_tool(args)
    
    args = [voyager_h5, '-S', '-p', 'ank', '-s', OUTDIR + 'test.png']
    cmd_tool(args)
    
    # Blank DC to .h5

    args = [voyager_h5, '-D', '-H', '-o', OUTDIR + 'test.h5']
    cmd_tool(args)

    # Blank DC to .fil

    args = [voyager_h5, '-D', '-F', '-o', OUTDIR + 'test.fil']
    cmd_tool(args)
    
    # info with foff negative

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
    print("\n===== test_cmd_arguments")
    args = [voyager_h5, '-H', '-F', '-o', OUTDIR + 'test.fil']
    with pytest.raises(ValueError):
        cmd_tool(args)

def test_neg_blank_dc():
    print("\n===== test_neg_blank_dc")
    wf = bl.Waterfall(voyager_h5)
    wf.blank_dc(0)
    wf.blank_dc(1.1)
    del wf

def test_get_chunk_dimensions():
    print("\n===== test_get_chunk_dimensions")
    wf = bl.Waterfall(voyager_h5)
    
    wf.header['foff'] = 0.99e-5
    assert wf._get_chunk_dimensions() == (1, 1, 1048576)
    wf.header['foff'] = 1.1e-5
    
    wf.header['tsamp'] = 0.99e-3
    assert wf._get_chunk_dimensions() == (2048, 1, 512)
    wf.header['tsamp'] = 1.1e-3
    
    wf.header['foff'] = 0.99e-2
    assert wf._get_chunk_dimensions() == (10, 1, 65536)
    
    wf.header['foff'] = 1e-1
    assert wf._get_chunk_dimensions() == (1, 1, 512)
    
    del wf

def test_neg_info_foff():
    print("\n===== test_neg_info_foff")
    wf = bl.Waterfall(voyager_h5)
    wf.header['foff'] = -1
    wf.info()
    del wf

def test_bug_no_filename():
    print("\n===== test_bug_no_filename")
    with pytest.raises(ValueError):
        bl.Waterfall()

