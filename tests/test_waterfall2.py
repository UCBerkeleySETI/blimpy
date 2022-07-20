import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
import blimpy as bl
from tests.data import voyager_h5

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
    wf_voya1 = bl.Waterfall(voyager_h5)
    wf_voya2 = bl.Waterfall(header_dict=wf_voya1.header, data_array=wf_voya1.data)
    wf_voya2.info()


if __name__ == "__main__":
    test_waterfall_stream_1()
    test_waterfall_stream_2()
