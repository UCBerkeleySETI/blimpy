from filterbank import Filterbank
import pylab as plt
import numpy as np

def test_voyager():
    filename = '/workdata/bl/data/voyager_f1032192_t300_v2.fil'
    
    fb = Filterbank(filename)
    fb.info()
    fb.plot_spectrum()
    plt.show()
    
    fb = Filterbank(filename, f_start=8420, f_stop=8420.5)
    fb.info()
    fb.plot_spectrum()
    plt.show()
    
def test_fix_header():
    filename = 'voyager_mod.fil'
    print read_header(filename)['ibeam']

    fix_header(filename, 'ibeam', 7)
    print read_header(filename)['ibeam']

    fix_header(filename, 'ibeam', 1)
    print read_header(filename)['ibeam']

    fix_header(filename, 'ibeam', 13)
    print read_header(filename)['ibeam']

    fix_header(filename, 'rawdatafile', './blc3_9bit_guppi_57386_VOYAGER1_0004.0000.raw')

if __name__ == "__main__":
    test_voyager()
