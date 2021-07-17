from blimpy import Waterfall
from tests.data import voyager_h5, test_ifs_fil


HIRES_THRESHOLD = 2**20


def test_ncc_chan_bw():
    wf = Waterfall(voyager_h5)
    print("test_bcc_chan_bw: telescope_id:", wf.header['telescope_id'])
    print("test_ncc_chan_bw: nchans:", wf.header['nchans'])

    n_coarse_chan = wf.calc_n_coarse_chan(16)
    print("test_ncc_chan_bw: n_coarse_chan [chan_bw=16]:", n_coarse_chan)
    assert n_coarse_chan == 64

    n_coarse_chan = wf.calc_n_coarse_chan(1)
    print("test_ncc_chan_bw: n_coarse_chan [chan_bw=1]:", n_coarse_chan)
    assert n_coarse_chan > 2.9 and n_coarse_chan < 3.0
    
    
def test_ncc_gbt():
    wf = Waterfall(test_ifs_fil)
    wf.header['telescope_id'] = 6
    print("test_ncc_gbt: telescope_id:", wf.header['telescope_id'])
    print("test_ncc_gbt: starting nchans:", wf.header['nchans'])

    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_gbt: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 64

    wf.header['nchans'] = 3 * HIRES_THRESHOLD
    print("\ntest_ncc_gbt: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_gbt: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 3

    wf.header['nchans'] = HIRES_THRESHOLD
    print("\ntest_ncc_gbt: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_gbt: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 1

    wf.header['nchans'] = HIRES_THRESHOLD - 1
    print("\ntest_ncc_gbt: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_gbt: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 64

    wf.header['nchans'] = HIRES_THRESHOLD + 1
    print("\ntest_ncc_gbt: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_gbt: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 64


def test_ncc_42():
    wf = Waterfall(test_ifs_fil)
    wf.header['telescope_id'] = 42
    print("test_ncc_42: telescope_id:", wf.header['telescope_id'])
    print("test_ncc_42: starting nchans:", wf.header['nchans'])

    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_42: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 64

    wf.header['nchans'] = 3 * HIRES_THRESHOLD
    print("\ntest_ncc_42: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_42: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 3

    wf.header['nchans'] = HIRES_THRESHOLD
    print("\ntest_ncc_42: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_42: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 1

    wf.header['nchans'] = HIRES_THRESHOLD - 1
    print("\ntest_ncc_42: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_42: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 64

    wf.header['nchans'] = HIRES_THRESHOLD + 1
    print("\ntest_ncc_42: nchans:", wf.header['nchans'])
    n_coarse_chan = wf.calc_n_coarse_chan()
    print("test_ncc_42: n_coarse_chan [chan_bw=None]:", n_coarse_chan)
    assert n_coarse_chan == 64


if __name__ == "__main__":
    wf = Waterfall(test_ifs_fil)
    wf.info()
    test_ncc_chan_bw()
    test_ncc_gbt()
    test_ncc_42()
