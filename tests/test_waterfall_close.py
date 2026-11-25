"""
Waterfall must release its input HDF5 file promptly (issue #177).

Before: the plotting methods were bound per instance in __init__
(six.create_bound_method(f, self)), so every Waterfall referenced itself.
Such an object - and the h5py.File it holds - is only freed by the cyclic
garbage collector, at an unpredictable later time. On Windows the file can't
be deleted, renamed or overwritten meanwhile (PermissionError, WinError 32).

Self-contained: the test file is written with blimpy itself (no downloads).
"""
import gc
import os

import h5py
import numpy as np
import pytest
from astropy.coordinates import Angle

import blimpy as bl


def _write_small_h5(path):
    nchans, ntints = 64, 8
    header = {"az_start": 0.0, "data_type": 1,
              "fch1": 8421.386717353016, "foff": -2.7939677238464355e-06,
              "ibeam": 1, "machine_id": 42, "nbeams": 1, "nbits": 32,
              "nchans": nchans, "nifs": 1, "rawdatafile": "nil",
              "source_name": "Test_Source",
              "src_raj": Angle("17:10:03.984 hours"),
              "src_dej": Angle("12:10:58.8 degrees"),
              "telescope_id": 42, "tstart": 57650.78209490741,
              "tsamp": 18.253611008, "za_start": 0.0}
    data = np.arange(ntints * nchans, dtype=np.float32).reshape(ntints, 1, nchans)
    bl.Waterfall(header_dict=header, data_array=data).write_to_hdf5(str(path))
    return path


def _open_handles(path):
    target = os.path.normcase(os.path.abspath(str(path)))
    n = 0
    for fid in h5py.h5f.get_obj_ids(types=h5py.h5f.OBJ_FILE):
        name = fid.name.decode() if isinstance(fid.name, bytes) else str(fid.name)
        if os.path.normcase(os.path.abspath(name)) == target:
            n += 1
    return n


@pytest.fixture
def small_h5(tmp_path):
    return _write_small_h5(tmp_path / "small.h5")


@pytest.fixture
def no_cyclic_gc():
    """Prove the file is released by plain reference counting, not by a lucky gc pass."""
    gc.collect()
    gc.disable()
    yield
    gc.enable()


@pytest.mark.parametrize("load_data", [False, True])
def test_dropping_a_waterfall_closes_its_file(small_h5, no_cyclic_gc, load_data):
    wf = bl.Waterfall(str(small_h5), load_data=load_data)
    assert _open_handles(small_h5) == 1
    del wf
    assert _open_handles(small_h5) == 0
    os.remove(small_h5)                       # raised WinError 32 on Windows before the fix
    assert not small_h5.exists()


def test_no_reference_cycle(small_h5):
    wf = bl.Waterfall(str(small_h5), load_data=False)
    assert not [k for k, v in vars(wf).items() if getattr(v, "__self__", None) is wf]
    wf.close()


def test_close_and_context_manager(small_h5):
    wf = bl.Waterfall(str(small_h5))
    data = wf.data.copy()
    wf.close()
    wf.close()                                # idempotent
    assert _open_handles(small_h5) == 0
    assert np.array_equal(wf.data, data)      # loaded data stays usable after close
    with bl.Waterfall(str(small_h5), load_data=False) as wf2:
        assert wf2.header["nchans"] == 64 and _open_handles(small_h5) == 1
    assert _open_handles(small_h5) == 0


def test_plot_methods_still_bound(small_h5):
    wf = bl.Waterfall(str(small_h5))
    for name in ("plot_spectrum", "plot_waterfall", "plot_kurtosis",
                 "plot_time_series", "plot_all", "plot_spectrum_min_max"):
        assert getattr(wf, name).__self__ is wf, name
    wf.close()
    # in-memory Waterfall too
    wf_mem = bl.Waterfall(header_dict=dict(wf.header), data_array=wf.data)
    assert wf_mem.plot_spectrum.__self__ is wf_mem
    wf_mem.close()                            # no file: a harmless no-op
