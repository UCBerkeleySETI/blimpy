from filterbank import Filterbank, read_header, fix_header
import pylab as plt
import numpy as np
import os
from pprint import pprint

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
    
def test_voyager_extract():
    filename = '/workdata/bl/data/voyager_f1032192_t300_v2.fil'
    new_filename = 'voyager_ext.fil'
    
    fb = Filterbank(filename, f_start=8420.1, f_stop=8420.3)
    fb.info()
    fb.plot_spectrum()
    plt.show()
    
    fb.write_to_filterbank(new_filename)
    
    fb2 = Filterbank(new_filename)
    fb2.info()
    fb2.plot_spectrum()  
    plt.show()
    
    os.remove(new_filename)  

def test_voyager_fix_header():
    filename = '/workdata/bl/data/voyager_f1032192_t300_v2.fil'
    new_filename = 'voyager_ext.fil'

    fb = Filterbank(filename, f_start=8420.1, f_stop=8420.3)
    fb.write_to_filterbank(new_filename)
    fb = Filterbank(new_filename)
    
    filename = new_filename
    assert read_header(filename)['ibeam'] == 1

    fix_header(filename, 'ibeam', 7)
    assert read_header(filename)['ibeam'] == 7

    fix_header(filename, 'ibeam', 1)
    assert read_header(filename)['ibeam'] == 1

    fix_header(filename, 'ibeam', 13)
    assert read_header(filename)['ibeam'] == 13
    
    pprint(read_header(filename))

    fix_header(filename, 'rawdatafile', './blc3_9bit_guppi_57386_VOYAGER1_0004.0000.raw')
    assert read_header(filename)['rawdatafile'] == './blc3_9bit_guppi_57386_VOYAGER1_0004.0000.raw'
    fix_header(filename, 'rawdatafile', './blc3_2bit_guppi_57386_VOYAGER1_0004.0000.raw')
    assert read_header(filename)['rawdatafile'] == './blc3_2bit_guppi_57386_VOYAGER1_0004.0000.raw'
    
    os.remove(new_filename)
    
if __name__ == "__main__":
    #test_voyager()
    #test_voyager_extract()
    test_voyager_fix_header()
