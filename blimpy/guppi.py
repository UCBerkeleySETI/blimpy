#!/usr/bin/env python
"""
# guppi.py

A python file handler for guppi RAW files from the GBT.

The guppi raw format consists of a FITS-like header, followed by a block of data,
and repeated over and over until the end of the file.
"""

import numpy as np
import os
import time
from pprint import pprint
from astropy.coordinates import Angle

from .utils import unpack, rebin
import sys

PYTHON3 = sys.version_info >= (3, 0)

# Check if $DISPLAY is set (for handling plotting on remote machines with no X-forwarding)
if 'DISPLAY' in os.environ.keys():

    try:
        import matplotlib
        #matplotlib.use('Qt5Agg')
    except ImportError:
        pass
    import pylab as plt
else:
    import matplotlib
    matplotlib.use('Agg')
    import pylab as plt

###
# Config values
###

MAX_PLT_POINTS      = 65536 * 4              # Max number of points in matplotlib plot
MAX_IMSHOW_POINTS   = (8192, 4096)           # Max number of points in imshow plot
MAX_DATA_ARRAY_SIZE = 1024 * 1024 * 1024     # Max size of data array to load into memory



class EndOfFileError(Exception):
    pass

class GuppiRaw(object):
    """ Python class for reading Guppi raw files

    Args:
        filename (str): name of the .raw file to open

    Optional args:
        n_blocks (int): if number of blocks to read is known, set it here.
                        This saves seeking through the file to check how many
                        integrations there are in the file.
    """
    def __init__(self, filename, n_blocks=None):
        self.filename = filename
        if PYTHON3:
            self.file_obj = open(filename, 'rb')
        else:
            self.file_obj = open(filename, 'r')
        self.filesize = os.path.getsize(filename)

        if not n_blocks:
            self.n_blocks = self.find_n_data_blocks()

        else:
            self.n_blocks = n_blocks

        self._d = np.zeros(1, dtype='complex64')
        self._d_x = np.zeros(1, dtype='int8')
        self._d_y = np.zeros(1, dtype='int8')

    def __enter__(self):
         return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.file_obj.close()

    def __repr__(self):
        return "<GuppiRaw file handler for %s>" % self.filename

    def read_header(self):
        """ Read next header (multiple headers in file)

        Returns:
            (header, data_idx) - a dictionary of keyword:value header data and
            also the byte index of where the corresponding data block resides.
        """
        start_idx = self.file_obj.tell()
        key, val = '', ''

        header_dict = {}
        keep_reading = True

        first_line = self.file_obj


        try:
            while keep_reading:
                if start_idx + 80 > self.filesize:
                    keep_reading = False
                    raise EndOfFileError
                line = self.file_obj.read(80)
                if PYTHON3:
                    line = line.decode("utf-8")
                #print line
                if line.startswith('END'):
                    keep_reading = False
                    break
                else:
                    key, val = line.split('=')
                    key, val = key.strip(), val.strip()

                    if "'" in val:
                        # Items in quotes are strings
                        val = str(val.strip("'").strip())
                    elif "." in val:
                        # Items with periods are floats (if not a string)
                        val = float(val)
                    else:
                        # Otherwise it's an integer
                        val = int(val)

                header_dict[key] = val
        except ValueError:
            print("CURRENT LINE: ", line)
            print("BLOCK START IDX: ",  start_idx)
            print("FILE SIZE: ",  self.filesize)
            print("NEXT 512 BYTES: \n")
            print(self.file_obj.read(512))
            raise

        data_idx = self.file_obj.tell()

        # Seek past padding if DIRECTIO is being used
        if "DIRECTIO" in header_dict.keys():
            if int(header_dict["DIRECTIO"]) == 1:
                if data_idx % 512:
                    data_idx += (512 - data_idx % 512)

        self.file_obj.seek(start_idx)
        return header_dict, data_idx

    def read_first_header(self):
        """ Read first header in file

        Returns:
            header (dict): keyword:value pairs of header metadata
        """
        self.file_obj.seek(0)
        header_dict, pos = self.read_header()
        self.file_obj.seek(0)
        return header_dict

    def read_next_data_block_shape(self):
        header, data_idx = self.read_header()
        n_chan = int(header['OBSNCHAN'])
        n_pol  = int(header['NPOL'])
        n_samples = int(header['BLOCSIZE']) / n_chan / n_pol
        
        is_chanmaj = False
        if 'CHANMAJ' in header.keys():
            if int(header['CHANMAJ']) == 1:
                is_chanmaj = True
        if is_chanmaj:
            dshape = (int(n_samples), n_chan, n_pol)
        else:
            dshape = (n_chan, int(n_samples), n_pol)
        return dshape

    def read_next_data_block_int8(self):
        """ Read the next block of data and its header

        Returns: (header, data)
            header (dict): dictionary of header metadata
            data (np.array): Numpy array of data, converted into to complex64.
        """
        header, data_idx = self.read_header()
        self.file_obj.seek(data_idx)

        # Read data and reshape

        n_chan = int(header['OBSNCHAN'])
        n_pol  = int(header['NPOL'])
        n_samples = int(header['BLOCSIZE']) / n_chan / n_pol
        n_bit = int(header['NBITS'])


        d = np.fromfile(self.file_obj, count=header['BLOCSIZE'], dtype='int8')

        # Handle 2-bit and 4-bit data
        if n_bit != 8:
            d = unpack(d, n_bit)

        d = d.reshape((n_chan, n_samples, n_pol))    # Real, imag

        if self._d_x.shape != d[..., 0:2].shape:
                self._d_x = np.ascontiguousarray(np.zeros(d[..., 0:2].shape, dtype='int8'))
                self._d_y = np.ascontiguousarray(np.zeros(d[..., 2:4].shape, dtype='int8'))

        self._d_x[:] = d[..., 0:2]
        self._d_y[:] = d[..., 2:4]
        return header, self._d_x, self._d_y
    

    def read_next_data_block_int8_2x(self):
        """ Read the next block of data and its header

        Returns: (header, data)
            header (dict): dictionary of header metadata
            data (np.array): Numpy array of data, converted into to complex64.
        """
        header, data_idx = self.read_header()
        self.file_obj.seek(data_idx)

        # Read data and reshape

        n_chan = int(header['OBSNCHAN'])
        n_pol  = int(header['NPOL'])
        n_samples = int(header['BLOCSIZE']) / n_chan / n_pol
        n_bit = int(header['NBITS'])

        d = np.fromfile(self.file_obj, count=header['BLOCSIZE'], dtype='int8')

        header, data_idx = self.read_header()
        self.file_obj.seek(data_idx)
        d2 = np.fromfile(self.file_obj, count=header['BLOCSIZE'], dtype='int8')

        # Handle 2-bit and 4-bit data
        if n_bit != 8:
            d = unpack(d, n_bit)

        d = d.reshape((n_chan, n_samples, n_pol))    # Real, imag
        d2 = d2.reshape((n_chan, n_samples, n_pol))
        d = np.concatenate((d, d2), axis=1)        
        print(d.shape)

        if self._d_x.shape != (n_chan, n_samples*2, n_pol):
                self._d_x = np.ascontiguousarray(np.zeros(d[..., 0:2].shape, dtype='int8'))
                self._d_y = np.ascontiguousarray(np.zeros(d[..., 2:4].shape, dtype='int8'))

        self._d_x[:] = d[..., 0:2]
        self._d_y[:] = d[..., 2:4]
        return header, self._d_x, self._d_y
    
    def read_next_data_block(self):
        """ Read the next block of data and its header

        Returns: (header, data)
            header (dict): dictionary of header metadata
            data (np.array): Numpy array of data, converted into to complex64.
        """
        header, data_idx = self.read_header()
        self.file_obj.seek(data_idx)

        # Read data and reshape

        n_chan = int(header['OBSNCHAN'])
        n_pol  = int(header['NPOL'])
        n_samples = int(header['BLOCSIZE']) / n_chan / n_pol
        n_bit = int(header['NBITS'])


        d = np.ascontiguousarray(np.fromfile(self.file_obj, count=header['BLOCSIZE'], dtype='int8'))

        # Handle 2-bit and 4-bit data
        if n_bit != 8:
            d = unpack(d, n_bit)
        
        dshape = self.read_next_data_block_shape()
        
        d = d.reshape(dshape)    # Real, imag

        if self._d.shape != d.shape:
            self._d = np.zeros(d.shape, dtype='float32')

        self._d[:] = d

        return header, self._d[:].view('complex64')
    

        

    def find_n_data_blocks(self):
        """ Seek through the file to find how many data blocks there are in the file

        Returns:
            n_blocks (int): number of data blocks in the file
        """
        self.file_obj.seek(0)
        header0, data_idx0 = self.read_header()

        self.file_obj.seek(data_idx0)
        block_size = int(header0['BLOCSIZE'])
        n_bits     = int(header0['NBITS'])
        self.file_obj.seek(int(header0['BLOCSIZE']), 1)
        n_blocks = 1
        end_found = False
        while not end_found:
            try:
                header, data_idx = self.read_header()
                self.file_obj.seek(data_idx)
                self.file_obj.seek(header['BLOCSIZE'], 1)
                n_blocks += 1
            except EndOfFileError:
                end_found = True
                break


        self.file_obj.seek(0)
        return n_blocks

    def reset_index(self):
        """ Return file_obj seek to start of file """
        self.file_obj.seek(0)

    def print_stats(self):
        """ Compute some basic stats on the next block of data """

        header, data = self.read_next_data_block()
        data = data.view('float32')

        print("AVG: %2.3f" % data.mean())
        print("STD: %2.3f" % data.std())
        print("MAX: %2.3f" % data.max())
        print("MIN: %2.3f" % data.min())

        import pylab as plt

    def plot_histogram(self, filename=None):
        """ Plot a histogram of data values """
        header, data = self.read_next_data_block()
        data = data.view('float32')

        plt.figure("Histogram")
        plt.hist(data.flatten(), 65, facecolor='#cc0000')
        if filename:
            plt.savefig(filename)
        plt.show()

    def plot_spectrum(self, filename=None, plot_db=True):
        """ Do a (slow) numpy FFT and take power of data """
        header, data = self.read_next_data_block()

        print("Computing FFT...")
        d_xx_fft = np.abs(np.fft.fft(data[..., 0]))
        d_xx_fft = d_xx_fft.flatten()

        # Rebin to max number of points
        dec_fac_x = 1
        if d_xx_fft.shape[0] > MAX_PLT_POINTS:
            dec_fac_x = d_xx_fft.shape[0] / MAX_PLT_POINTS

        d_xx_fft  = rebin(d_xx_fft, dec_fac_x)

        print("Plotting...")
        if plot_db:
            plt.plot(10 * np.log10(d_xx_fft))
            plt.ylabel("Power [dB]")
        else:
            plt.plot(d_xx_fft)
            plt.ylabel("Power")
            plt.xlabel("Channel")
        plt.title(self.filename)
        if filename:
            plt.savefig(filename)
        plt.show()

    def generate_filterbank_header(self, nchans=1, ):
        """ Generate a blimpy header dictionary """
        gp_head = self.read_first_header()
        fb_head = {}

        telescope_str = gp_head.get("TELESCOP", "unknown")
        if telescope_str in ('GBT', 'GREENBANK'):
            fb_head["telescope_id"] = 6
        elif telescope_str in ('PKS', 'PARKES'):
            fb_head["telescop_id"] = 7
        else:
            fb_head["telescop_id"] = 0

        # Using .get() method allows us to fill in default values if not present
        fb_head["source_name"] = gp_head.get("SRC_NAME", "unknown")
        fb_head["az_start"]    = gp_head.get("AZ", 0)
        fb_head["za_start"]    = gp_head.get("ZA", 0)

        fb_head["src_raj"]     = Angle(str(gp_head.get("RA", 0.0)) + "hr")
        fb_head["src_dej"]     = Angle(str(gp_head.get("DEC", 0.0)) + "deg")
        fb_head["rawdatafile"] = self.filename

        # hardcoded
        fb_head["machine_id"]    = 20
        fb_head["data_type"]     = 1      # blio datatype
        fb_head["barycentric"]      = 0
        fb_head["pulsarcentric"] = 0
        fb_head["nbits"]         = 32

        # TODO - compute these values. Need to figure out the correct calcs
        fb_head["tstart"]        = 0.0
        fb_head["tsamp"]         = 1.0
        fb_head["fch1"]          = 0.0
        fb_head["foff"]          = 187.5 / nchans

        # Need to be updated based on output specs
        fb_head["nchans"]        = nchans
        fb_head["nifs"]          = 1
        fb_head["nbeams"]        = 1

        return fb_head


def cmd_tool(args=None):
    """ Command line tool for plotting and viewing info on guppi raw files """
    
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Command line utility for creating spectra from GuppiRaw files.")

    parser.add_argument('filename', type=str, help='Name of file to read')
    parser.add_argument('-o', dest='outdir', type=str, default='./', help='output directory for PNG files')
    args = parser.parse_args()

    r = GuppiRaw(args.filename)

    r.print_stats()
    bname = os.path.splitext(os.path.basename(args.filename))[0]
    bname = os.path.join(args.outdir, bname)
    r.plot_histogram(filename="%s_hist.png" % bname)
    r.plot_spectrum(filename="%s_spec.png" % bname)


if __name__ == "__main__":
    cmd_tool()
