import os
import numpy as np
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
import blimpy as bl
from blimpy.plotting import plot_waterfall, plot_spectrum_min_max
from tests.data import voyager_h5, test_h5, test_fil


TEST_DATA_DIR = os.path.dirname(voyager_h5)
RTOL = 1e-05


def compare_hdr_fields(hdr1, hdr2, fieldname):
    if isinstance(hdr1[fieldname], float):
        if np.isclose(hdr1[fieldname], hdr2[fieldname], rtol=RTOL):
            return 0
    else: # not a float: int or str
        if hdr1[fieldname] == hdr2[fieldname]:
            return 0
    print(f"*** compare_hdr_fields: {hdr1[fieldname]} != {hdr2[fieldname]}")
    return 1


def compare_data_vectors(label, vec1, vec2):
    result = np.isclose(vec1, vec2, rtol=RTOL)
    if False in result:
        print(f"*** compare_data_vectors: {label}: {vec1} != {vec2}")
        return 1
    return 0


def spot_check_data(data1, data2, n_ints_in_file, n_channels_in_file):
    nw1 = data1[0, 0, 0:3]
    ne1 = data1[0, 0, -4:-1]
    sw1 = data1[n_ints_in_file - 1, 0, 0:3]
    se1 = data1[n_ints_in_file - 1, 0, -4:-1]
    centre_row = n_ints_in_file // 2
    centre_col = n_channels_in_file // 2
    bullseye1 = data1[centre_row, 0, centre_col - 1 : centre_col + 2]

    nw2 = data2[0, 0, 0:3]
    ne2 = data2[0, 0, -4:-1]
    sw2 = data2[n_ints_in_file - 1, 0, 0:3]
    se2 = data2[n_ints_in_file - 1, 0, -4:-1]
    centre_row = n_ints_in_file // 2
    centre_col = n_channels_in_file // 2
    bullseye2 = data2[centre_row, 0, centre_col - 1 : centre_col + 2]

    n_errors = 0
    n_errors += compare_data_vectors("nw", nw1, nw2)
    n_errors += compare_data_vectors("ne", ne1, ne2)
    n_errors += compare_data_vectors("sw", sw1, sw2)
    n_errors += compare_data_vectors("se", se1, se2)
    n_errors += compare_data_vectors("bullseye", bullseye1, bullseye2)
    return n_errors


def compare_headers(hdr1, hdr2):
    n_errors = 0
    n_errors += compare_hdr_fields(hdr1, hdr2, "fch1")
    n_errors += compare_hdr_fields(hdr1, hdr2, "nchans")
    n_errors += compare_hdr_fields(hdr1, hdr2, "nifs")
    n_errors += compare_hdr_fields(hdr1, hdr2, "nbits")
    n_errors += compare_hdr_fields(hdr1, hdr2, "source_name")
    n_errors += compare_hdr_fields(hdr1, hdr2, "telescope_id")
    n_errors += compare_hdr_fields(hdr1, hdr2, "tsamp")
    n_errors += compare_hdr_fields(hdr1, hdr2, "tstart")
    n_errors += compare_hdr_fields(hdr1, hdr2, "src_raj")
    n_errors += compare_hdr_fields(hdr1, hdr2, "src_dej")
    return n_errors


def test_waterfall_stream_1():

    print("\n===== test_waterfall_stream_1")

    source_name = "Not_Voyager_1"
    src_raj = Angle("17:10:03.984 hours")
    src_dej = Angle("12:10:58.8 degrees")
    tstart = 57650.78209490741
    tsamp = 18.253611008
    f_start = 8418.457032646984
    f_stop = 8421.386717353016
    n_fine_chans = 20
    n_tints = 8

    foff = (f_stop - f_start) / float(n_fine_chans)

    header = {"az_start": 0.0, "data_type": 1,
                   "fch1": f_start, "foff": foff,
                   "ibeam": 1, "machine_id": 42, "nbeams": 1, "nbits": 32,
                   "nchans": n_fine_chans, "nifs": 1, "rawdatafile": "nil",
                   "source_name": source_name, "src_raj": src_raj, "src_dej": src_dej,
                   "telescope_id": 42,
                   "tstart": tstart, "tsamp": tsamp, "zs_tart": 0.0}

    data_matrix = np.zeros((n_tints, 1, n_fine_chans), dtype=np.float32)

    wf = bl.Waterfall(header_dict=header, data_array=data_matrix)
    print("\nwf:", wf)
    wf.info()


def test_waterfall_stream_2():

    print("\n===== test_waterfall_stream_2")
    wf_storage = bl.Waterfall(voyager_h5)
    wf_stream = bl.Waterfall(header_dict=wf_storage.header, data_array=wf_storage.data)
    assert compare_headers(wf_storage.header, wf_stream.header) == 0
    assert spot_check_data(wf_storage.data,
                           wf_stream.data,
                           wf_storage.n_ints_in_file,
                           wf_storage.n_channels_in_file) == 0
    wf_stream.info()
    plt.figure("Voyager 1", figsize=(10, 3))
    plt.subplot(1, 2, 1) # nrows=3, ncols=2, index=1 relative to 1
    plot_waterfall(wf_stream)
    plt.subplot(1, 2, 2) # nrows=3, ncols=2, index=2 relative to 1
    plot_spectrum_min_max(wf_stream)
    plt.tight_layout()
    plt.savefig(TEST_DATA_DIR + "/test_waterfall_stream_2.png")
    wf_stream.write_to_hdf5(test_h5)
    wf_stream.write_to_fil(test_fil)


if __name__ == "__main__":
    test_waterfall_stream_1()
    test_waterfall_stream_2()
