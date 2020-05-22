import numpy as np

import blimpy as bl
from tests.data import voyager_raw, voyager_block1

def test_guppi():
    gr = bl.guppi.GuppiRaw(voyager_raw)
    h1, data_block_x1, data_block_y1 = gr.read_next_data_block_int8()
    h2, data_block_x2, data_block_y2 = gr.read_next_data_block_int8()

    assert not np.array_equal(data_block_x1, data_block_x2) and not np.array_equal(data_block_y1, data_block_y2) \
        , "Data read from two blocks should not be equal"

    gr = bl.guppi.GuppiRaw(voyager_raw)
    h1, data_block_1 = gr.read_next_data_block()

    data_block_reference_1 = np.load(voyager_block1)

    assert np.array_equal(data_block_1[:, :1000, :], data_block_reference_1) \
        , "Data read should be consistent with previous versions"

    data_block_casted_1 = np.append(data_block_x1, data_block_y1, axis=2).astype('float32').view('complex64')

    assert np.array_equal(data_block_1, data_block_casted_1) \
        , "Reading as int8 then casting should be equal to reading directly as complex64"

# We have to keep instantiating objects because
# the plotting routines read data in a manner
# destructive to the object.    

def test_spectrum():
    gr = bl.guppi.GuppiRaw(voyager_raw)
    gr.plot_spectrum()

def test_histogram():
    gr = bl.guppi.GuppiRaw(voyager_raw)
    gr.plot_histogram()

def test_statistics():
    gr = bl.guppi.GuppiRaw(voyager_raw)
    gr.print_stats()

def test_fil_header():
    gr = bl.guppi.GuppiRaw(voyager_raw)
    gr.generate_filterbank_header()

if __name__ == "__main__":
    test_guppi()
    print("OK")
