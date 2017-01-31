#!/usr/bin/env python
"""
# filterbank.py

Python class and command line utility for reading and plotting filterbank files.

This provides a class, Filterbank(), which can be used to read a .fil file:

    ````
    fil = Filterbank('test_psr.fil')
    print fil.header
    print fil.data.shape
    print fil.freqs

    plt.figure()
    fil.plot_spectrum(t=0)
    plt.show()
    ````

TODO: check the file seek logic works correctly for multiple IFs

"""

import os
import sys
import struct
import numpy as np
from pprint import pprint

from astropy import units as u
from astropy.coordinates import Angle
import scipy.stats
from matplotlib.ticker import NullFormatter

from utils import db, lin, rebin, closest

try:
    import h5py
    HAS_HDF5 = True
except ImportError:
    HAS_HDF5 = False

# Check if $DISPLAY is set (for handling plotting on remote machines with no X-forwarding)
if os.environ.has_key('DISPLAY'):
    import pylab as plt
else:
    import matplotlib
    matplotlib.use('Agg')
    import pylab as plt

###
# Config values
###

MAX_PLT_POINTS      = 65536                  # Max number of points in matplotlib plot
MAX_IMSHOW_POINTS   = (8192, 4096)           # Max number of points in imshow plot
MAX_DATA_ARRAY_SIZE = 1024 * 1024 * 1024     # Max size of data array to load into memory
MAX_HEADER_BLOCKS   = 100                    # Max size of header (in 512-byte blocks)

###
# Header parsing
###

# Dictionary of allowed keywords and their types
# Here are the keywordss that a filter bank file may
# contain.  Items marked with "[*]" are not yet # supported.  See docs for
# indivisuabl attribtues for more detailed info.
#
#   * telescope_id (int): 0=fake data; 1=Arecibo; 2=Ooty... others to be added
#   * machine_id (int): 0=FAKE; 1=PSPM; 2=WAPP; 3=OOTY... others to be added
#   * data_type (int): 1=filterbank; 2=time series... others to be added
#   * rawdatafile (string): the name of the original data file
#   * source_name (string): the name of the source being observed by the telescope
#   * barycentric (int): equals 1 if data are barycentric or 0 otherwise
#   * pulsarcentric (int): equals 1 if data are pulsarcentric or 0 otherwise
#   * az_start (double): telescope azimuth at start of scan (degrees)
#   * za_start (double): telescope zenith angle at start of scan (degrees)
#   * src_raj (double): right ascension (J2000) of source (hours, converted from hhmmss.s)
#   * src_dej (double): declination (J2000) of source (degrees, converted from ddmmss.s)
#   * tstart (double): time stamp (MJD) of first sample
#   * tsamp (double): time interval between samples (s)
#   * nbits (int): number of bits per time sample
#   * nsamples (int): number of time samples in the data file (rarely used any more)
#   * fch1 (double): centre frequency (MHz) of first filterbank channel
#   * foff (double): filterbank channel bandwidth (MHz)
#   * FREQUENCY_START [*] (character): start of frequency table (see below for explanation)
#   * fchannel [*] (double): frequency channel value (MHz)
#   * FREQUENCY_END [*] (character): end of frequency table (see below for explanation)
#   * nchans (int): number of filterbank channels
#   * nifs (int): number of seperate IF channels
#   * refdm (double): reference dispersion measure (pc/cm**3)
#   * period (double): folding period (s)
#   * nbeams (int):total number of beams (?)
#   * ibeam (int): number of the beam in this file (?)

header_keyword_types = {
    'telescope_id' : '<l',
    'machine_id'   : '<l',
    'data_type'    : '<l',
    'barycentric'  : '<l',
    'pulsarcentric': '<l',
    'nbits'        : '<l',
    'nsamples'     : '<l',
    'nchans'       : '<l',
    'nifs'         : '<l',
    'nbeams'       : '<l',
    'ibeam'        : '<l',
    'rawdatafile'  : 'str',
    'source_name'  : 'str',
    'az_start'     : '<d',
    'za_start'     : '<d',
    'tstart'       : '<d',
    'tsamp'        : '<d',
    'fch1'         : '<d',
    'foff'         : '<d',
    'refdm'        : '<d',
    'period'       : '<d',
    'src_raj'      : 'angle',
    'src_dej'      : 'angle',
    }

def grab_header(filename):
    """ Extract the filterbank header from the file

    Args:
        filename (str): name of file to open

    Returns:
        header_str (str): filterbank header as a binary string
    """
    f = open(filename, 'rb')
    eoh_found = False

    header_str = ''
    header_sub_count = 0
    while not eoh_found:
        header_sub = f.read(512)
        header_sub_count += 1
        if 'HEADER_START' in header_sub:
            idx_start = header_sub.index('HEADER_START') + len('HEADER_START')
            header_sub = header_sub[idx_start:]

        if 'HEADER_END' in header_sub:
            eoh_found = True
            idx_end = header_sub.index('HEADER_END')
            header_sub = header_sub[:idx_end]

        if header_sub_count >= MAX_HEADER_BLOCKS:
            raise RuntimeError("MAX HEADER LENGTH REACHED. THIS FILE IS FUBARRED.")
        header_str += header_sub

    f.close()
    return header_str

def len_header(filename):
    """ Return the length of the filterbank header, in bytes

    Args:
        filename (str): name of file to open

    Returns:
        idx_end (int): length of header, in bytes
    """
    with  open(filename, 'rb') as f:
        header_sub_count = 0
        eoh_found = False
        while not eoh_found:
            header_sub = f.read(512)
            header_sub_count += 1
            if 'HEADER_END' in header_sub:
                idx_end = header_sub.index('HEADER_END') + len('HEADER_END')
                eoh_found = True
                break

        idx_end = (header_sub_count -1) * 512 + idx_end
    return idx_end

def parse_header(filename):
    """ Parse a header of a filterbank, looking for allowed keywords

    Uses header_keyword_types dictionary as a lookup for data types.

    Args:
        filename (str): name of file to open

    Returns:
        header_dict (dict): A dictioary of header key:value pairs
    """
    header = grab_header(filename)
    header_dict = {}

    #print header
    for keyword in header_keyword_types.keys():
        if keyword in header:
            dtype = header_keyword_types.get(keyword, 'str')
            idx = header.index(keyword) + len(keyword)
            dtype = header_keyword_types[keyword]
            if dtype == '<l':
                val = struct.unpack(dtype, header[idx:idx+4])[0]
                header_dict[keyword] = val
            if dtype == '<d':
                val = struct.unpack(dtype, header[idx:idx+8])[0]
                header_dict[keyword] = val
            if dtype == 'str':
                str_len = struct.unpack('<L', header[idx:idx+4])[0]
                str_val = header[idx+4:idx+4+str_len]
                header_dict[keyword] = str_val
            if dtype == 'angle':
                val = struct.unpack('<d', header[idx:idx+8])[0]
                val = fil_double_to_angle(val)

                if keyword == 'src_raj':
                    val = Angle(val, unit=u.hour)
                else:
                    val = Angle(val, unit=u.deg)
                header_dict[keyword] = val

    return header_dict

def read_next_header_keyword(fh):
    """

    Args:
        fh (file): file handler

    Returns:
    """
    n_bytes = np.fromstring(fh.read(4), dtype='uint32')[0]
    #print n_bytes

    if n_bytes > 255:
        n_bytes = 16

    keyword = fh.read(n_bytes)

    #print keyword

    if keyword == 'HEADER_START' or keyword == 'HEADER_END':
        return keyword, 0, fh.tell()
    else:
        dtype = header_keyword_types[keyword]
        #print dtype
        idx = fh.tell()
        if dtype == '<l':
            val = struct.unpack(dtype, fh.read(4))[0]
        if dtype == '<d':
            val = struct.unpack(dtype, fh.read(8))[0]
        if dtype == 'str':
            str_len = np.fromstring(fh.read(4), dtype='int32')[0]
            val = fh.read(str_len)
        if dtype == 'angle':
            val = struct.unpack('<d', fh.read(8))[0]
            val = fil_double_to_angle(val)
            if keyword == 'src_raj':
                val = Angle(val, unit=u.hour)
            else:
                val = Angle(val, unit=u.deg)
        return keyword, val, idx

def read_header(filename, return_idxs=False):
    """ Read filterbank header and return a Python dictionary of key:value pairs

    Args:
        filename (str): name of file to open

    Optional args:
        return_idxs (bool): Default False. If true, returns the file offset indexes
                            for values

    returns

    """
    with open(filename, 'rb') as fh:
        header_dict = {}
        header_idxs = {}

        # Check this is a filterbank file
        keyword, value, idx = read_next_header_keyword(fh)

        try:
            assert keyword == 'HEADER_START'
        except AssertionError:
            raise RuntimeError("Not a valid filterbank file.")

        while True:
            keyword, value, idx = read_next_header_keyword(fh)
            if keyword == 'HEADER_END':
                break
            else:
                header_dict[keyword] = value
                header_idxs[keyword] = idx

    if return_idxs:
        return header_idxs
    else:
        return header_dict

def fix_header(filename, keyword, new_value):
    """ Apply a quick patch-up to a Filterbank header by overwriting a header value


    Args:
        filename (str): name of file to open and fix. WILL BE MODIFIED.
        keyword (stt):  header keyword to update
        new_value (long, double, angle or string): New value to write.

    Notes:
        This will overwrite the current value of the filterbank with a desired
        'fixed' version. Note that this has limited support for patching
        string-type values - if the length of the string changes, all hell will
        break loose.

    """

    # Read header data and return indexes of data offsets in file
    hd = read_header(filename)
    hi = read_header(filename, return_idxs=True)
    idx = hi[keyword]

    # Find out the datatype for the given keyword
    dtype = header_keyword_types[keyword]
    dtype_to_type = {'<l'  : np.int32,
                     'str' : str,
                     '<d'  : np.float64,
                     'angle' : to_sigproc_angle}
    value_dtype = dtype_to_type[dtype]

    # Generate the new string
    if value_dtype is str:
        if len(hd[keyword]) == len(new_value):
            val_str = np.int32(len(new_value)).tostring() + new_value
        else:
            raise RuntimeError("String size mismatch. Cannot update without rewriting entire file.")
    else:
        val_str = value_dtype(new_value).tostring()

    # Write the new string to file
    with open(filename, 'rb+') as fh:
        fh.seek(idx)
        fh.write(val_str)

def fil_double_to_angle(angle):
      """ Reads a little-endian double in ddmmss.s (or hhmmss.s) format and then
      converts to Float degrees (or hours).  This is primarily used to read
      src_raj and src_dej header values. """

      negative = (angle < 0.0)
      angle = np.abs(angle)

      dd = np.floor((angle / 10000))
      angle -= 10000 * dd
      mm = np.floor((angle / 100))
      ss = angle - 100 * mm
      dd += mm/60.0 + ss/3600.0

      if negative:
          dd *= -1

      return dd

###
# sigproc writing functions
###

def to_sigproc_keyword(keyword, value=None):
    """ Generate a serialized string for a sigproc keyword:value pair

    If value=None, just the keyword will be written with no payload.
    Data type is inferred by keyword name (via a lookup table)

    Args:
        keyword (str): Keyword to write
        value (None, float, str, double or angle): value to write to file

    Returns:
        value_str (str): serialized string to write to file.
    """
    if not value:
        return np.int32(len(keyword)).tostring() + keyword
    else:
        dtype = header_keyword_types[keyword]

        dtype_to_type = {'<l'  : np.int32,
                         'str' : str,
                         '<d'  : np.float64,
                         'angle' : to_sigproc_angle}

        value_dtype = dtype_to_type[dtype]

        if value_dtype is str:
            return np.int32(len(keyword)).tostring() + keyword + np.int32(len(value)).tostring() + value
        else:
            return np.int32(len(keyword)).tostring() + keyword + value_dtype(value).tostring()

def generate_sigproc_header(f):
    """ Generate a serialzed sigproc header which can be written to disk.

    Args:
        f (Filterbank object): Filterbank object for which to generate header

    Returns:
        header_str (str): Serialized string corresponding to header
    """

    header_string = ''
    header_string += to_sigproc_keyword('HEADER_START')

    for keyword in f.header.keys():
            if keyword == 'src_raj':
                header_string += to_sigproc_keyword('src_raj')  + to_sigproc_angle(f.header['src_raj'])
            elif keyword == 'src_dej':
                header_string += to_sigproc_keyword('src_dej')  + to_sigproc_angle(f.header['src_dej'])
            elif keyword == 'az_start' or keyword == 'za_start':
                header_string += to_sigproc_keyword(keyword)  + np.float64(f.header[keyword]).tostring()
            else:
                header_string += to_sigproc_keyword(keyword, f.header[keyword])

    header_string += to_sigproc_keyword('HEADER_END')
    return header_string

def to_sigproc_angle(angle_val):
    """ Convert an astropy.Angle to the ridiculous sigproc angle format string. """
    x = str(angle_val)

    if 'h' in x:
        d, m, s, ss = int(x[0:x.index('h')]), int(x[x.index('h')+1:x.index('m')]), \
        int(x[x.index('m')+1:x.index('.')]), float(x[x.index('.'):x.index('s')])
    if 'd' in x:
        d, m, s, ss = int(x[0:x.index('d')]), int(x[x.index('d')+1:x.index('m')]), \
        int(x[x.index('m')+1:x.index('.')]), float(x[x.index('.'):x.index('s')])
    num = str(d).zfill(2) + str(m).zfill(2) + str(s).zfill(2)+ '.' + str(ss).split(".")[-1]
    return np.float64(num).tostring()

###
# Main filterbank class
###

class Filterbank(object):
    """ Class for loading and plotting filterbank data """

    def __repr__(self):
        return "Filterbank data: %s" % self.filename

    def __init__(self, filename=None, f_start=None, f_stop=None,
                 t_start=None, t_stop=None, load_data=True,
                 header_dict=None, data_array=None):
        """ Class for loading and plotting filterbank data.

        This class parses the filterbank file and stores the header and data
        as objects:
            fb = Filterbank('filename_here.fil')
            fb.header        # filterbank header, as a dictionary
            fb.data          # filterbank data, as a numpy array

        Args:
            filename (str): filename of filterbank file.
            f_start (float): start frequency in MHz
            f_stop (float): stop frequency in MHz
            t_start (int): start integration ID
            t_stop (int): stop integration ID
            load_data (bool): load data. If set to False, only header will be read.
            header_dict (dict): Create filterbank from header dictionary + data array
            data_array (np.array): Create filterbank from header dict + data array
        """

        if filename:
            self.filename = filename
            if HAS_HDF5:
                if h5py.is_hdf5(filename):
                    self.read_hdf5(filename, f_start, f_stop, t_start, t_stop, load_data)
                else:
                    self.read_filterbank(filename, f_start, f_stop, t_start, t_stop, load_data)
            else:
                self.read_filterbank(filename, f_start, f_stop, t_start, t_stop, load_data)
        elif header_dict is not None and data_array is not None:
            self.gen_from_header(header_dict, data_array)
        else:
            pass

    def gen_from_header(self, header_dict, data_array, f_start=None, f_stop=None,
                        t_start=None, t_stop=None, load_data=True):
        self.filename = ''
        self.header = header_dict
        self.data = data_array
        self.n_ints_in_file = 0

        self._setup_freqs()



    def read_hdf5(self, filename, f_start=None, f_stop=None,
                        t_start=None, t_stop=None, load_data=True):
        self.header = {}
        self.filename = filename
        self.h5 = h5py.File(filename)
        for key, val in self.h5['data'].attrs.items():
            if key == 'src_raj':
                self.header[key] = Angle(val, unit='hr')
            elif key == 'src_dej':
                self.header[key] = Angle(val, unit='deg')
            else:
                self.header[key] = val

        self.data = self.h5["data"][:]
        self._setup_freqs()

        self.n_ints_in_file  = self.data.shape[0]
        self.file_size_bytes = os.path.getsize(self.filename)

    def _setup_freqs(self, f_start=None, f_stop=None):
        ## Setup frequency axis
        f0 = self.header['fch1']
        f_delt = self.header['foff']

        i_start, i_stop = 0, self.header['nchans']
        if f_start:
            i_start = (f_start - f0) / f_delt
        if f_stop:
            i_stop  = (f_stop - f0)  / f_delt

        #calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx  = np.int(i_stop)

        #create freq array
        if i_start < i_stop:
            i_vals = np.arange(chan_start_idx, chan_stop_idx)
        else:
            i_vals = np.arange(chan_stop_idx, chan_start_idx)

        self.freqs = f_delt * i_vals + f0

        if f_delt < 0:
            self.freqs = self.freqs[::-1]

        return i_start, i_stop, chan_start_idx, chan_stop_idx


    def read_filterbank(self, filename=None, f_start=None, f_stop=None,
                        t_start=None, t_stop=None, load_data=True):

        if filename is None:
            filename = self.filename

        self.header = read_header(filename)

        ## Setup frequency axis
        f0 = self.header['fch1']
        f_delt = self.header['foff']

        # keep this seperate!
        # file_freq_mapping =  np.arange(0, self.header['nchans'], 1, dtype='float64') * f_delt + f0

        #convert input frequencies into what their corresponding index would be

        i_start, i_stop, chan_start_idx, chan_stop_idx = self._setup_freqs(f_start, f_stop)

        n_bytes  = self.header['nbits'] / 8
        n_chans = self.header['nchans']
        n_chans_selected = self.freqs.shape[0]
        n_ifs   = self.header['nifs']


        # Load binary data
        self.idx_data = len_header(filename)
        f = open(filename, 'rb')
        f.seek(self.idx_data)
        filesize = os.path.getsize(self.filename)
        n_bytes_data = filesize - self.idx_data
        n_ints_in_file = n_bytes_data / (n_bytes * n_chans * n_ifs)

        # now check to see how many integrations requested
        ii_start, ii_stop = 0, n_ints_in_file
        if t_start:
            ii_start = t_start
        if t_stop:
            ii_stop = t_stop
        n_ints = ii_stop - ii_start

        # Seek to first integration
        f.seek(ii_start * n_bytes * n_ifs * n_chans, 1)

        # Set up indexes used in file read (taken out of loop for speed)
        i0 = np.min((chan_start_idx, chan_stop_idx))
        i1 = np.max((chan_start_idx, chan_stop_idx))

        if load_data:

            if n_ints * n_ifs * n_chans_selected > MAX_DATA_ARRAY_SIZE:
                print "Error: data array is too large to load. Either select fewer"
                print "points or manually increase MAX_DATA_ARRAY_SIZE."
                exit()

            self.data = np.zeros((n_ints, n_ifs, n_chans_selected), dtype='float32')

            for ii in range(n_ints):
                """d = f.read(n_bytes * n_chans * n_ifs)
                """

                for jj in range(n_ifs):

                    f.seek(n_bytes * i0, 1) # 1 = from current location
                    #d = f.read(n_bytes * n_chans_selected)
                    #bytes_to_read = n_bytes * n_chans_selected
                    if n_bytes == 4:
                        dd = np.fromfile(f, count=n_chans_selected, dtype='float32')
                    elif n_bytes == 2:
                        dd = np.fromfile(f, count=n_chans_selected, dtype='int16')
                    elif n_bytes == 1:
                        dd = np.fromfile(f, count=n_chans_selected, dtype='int8')

                    # Reverse array if frequency axis is flipped
                    if f_delt < 0:
                        dd = dd[::-1]

                    self.data[ii, jj] = dd

                    f.seek(n_bytes * (n_chans - i1), 1)  # Seek to start of next block
        else:
            print "Skipping data load..."
            self.data = np.array([0])

        # Finally add some other info to the class as objects
        self.n_ints_in_file  = n_ints_in_file
        self.file_size_bytes = filesize

        ## Setup time axis
        t0 = self.header['tstart']
        t_delt = self.header['tsamp']
        self.timestamps = np.arange(0, n_ints) * t_delt / 24./60./60 + t0

    def blank_dc(self, n_coarse_chan):
        """ Blank DC bins in coarse channels.

        Note: currently only works if entire filterbank file is read
        """
        n_chan = self.data.shape[2]
        n_chan_per_coarse = n_chan / n_coarse_chan

        mid_chan = n_chan_per_coarse / 2

        for ii in range(0, n_coarse_chan-1):
            ss = ii*n_chan_per_coarse
            self.data[..., ss+mid_chan-1] = self.data[..., ss+mid_chan]


    def info(self):
        """ Print header information """

        for key, val in self.header.items():
            if key == 'src_raj':
                val = val.to_string(unit=u.hour, sep=':')
            if key == 'src_dej':
                val = val.to_string(unit=u.deg, sep=':')
            print "%16s : %32s" % (key, val)

        print "\n%16s : %32s" % ("Num ints in file", self.n_ints_in_file)
        print "%16s : %32s" % ("Data shape", self.data.shape)
        print "%16s : %32s" % ("Start freq (MHz)", self.freqs[0])
        print "%16s : %32s" % ("Stop freq (MHz)", self.freqs[-1])

    def generate_freqs(self, f_start, f_stop):
        """
        returns frequency array [f_start...f_stop]
        """

        fch1 = self.header['fch1']
        foff = self.header['foff']

        #convert input frequencies into what their corresponding index would be
        i_start = (f_start - fch1) / foff
        i_stop  = (f_stop - fch1)  / foff

        #calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx  = np.int(i_stop)

        #create freq array
        i_vals = np.arange(chan_stop_idx, chan_start_idx, 1)

        freqs = foff * i_vals + fch1
        return freqs[::-1]

    def grab_data(self, f_start=None, f_stop=None, if_id=0):
        """ Extract a portion of data by frequency range.

        Args:
            f_start (float): start frequency in MHz
            f_stop (float): stop frequency in MHz
            if_id (int): IF input identification (req. when multiple IFs in file)

        Returns:
            (freqs, data) (np.arrays): frequency axis in MHz and data subset
        """
        i_start, i_stop = 0, None

        if f_start:
            i_start = closest(self.freqs, f_start)
        if f_stop:
            i_stop = closest(self.freqs, f_stop)

        plot_f    = self.freqs[i_start:i_stop]
        plot_data = self.data[:, if_id, i_start:i_stop]
        return plot_f, plot_data

    def plot_spectrum(self, t=0, f_start=None, f_stop=None, logged=False, if_id=0, c=None, **kwargs):
        """ Plot frequency spectrum of a given file

        Args:
            t (int): integration number to plot (0 -> len(data))
            logged (bool): Plot in linear (False) or dB units (True)
            if_id (int): IF identification (if multiple IF signals in file)
            c: color for line
            kwargs: keyword args to be passed to matplotlib plot()
        """
        ax = plt.gca()

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        if isinstance(t, int):
            print "extracting integration %i..." % t
            plot_data = plot_data[t]
        elif t == 'all':
            print "averaging along time axis..."
            plot_data = plot_data.mean(axis=0)
        else:
            raise RuntimeError("Unknown integration %s" % t)

        # Rebin to max number of points
        dec_fac_x = 1
        if plot_data.shape[0] > MAX_PLT_POINTS:
            dec_fac_x = plot_data.shape[0] / MAX_PLT_POINTS

        plot_data = rebin(plot_data, dec_fac_x, 1)
        plot_f    = rebin(plot_f, dec_fac_x, 1)

        if not c:
            kwargs['c'] = '#333333'

        if logged:
            plt.plot(plot_f, db(plot_data),label='Stokes I', **kwargs)
            plt.ylabel("Power [dB]")
        else:

            plt.plot(plot_f, plot_data,label='Stokes I', **kwargs)
            plt.ylabel("Power [counts]")
        plt.xlabel("Frequency [MHz]")
        plt.legend()

        try:
            plt.title(self.header['source_name'])
        except KeyError:
            plt.title(self.filename)

        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.xlim(plot_f[0], plot_f[-1])

    def plot_spectrum_min_max(self, t=0, f_start=None, f_stop=None, logged=False, if_id=0, c=None, **kwargs):
        """ Plot frequency spectrum of a given file

        Args:
            logged (bool): Plot in linear (False) or dB units (True)
            if_id (int): IF identification (if multiple IF signals in file)
            c: color for line
            kwargs: keyword args to be passed to matplotlib plot()
        """
        ax = plt.gca()

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        fig_max = plot_data[0].max()
        fig_min = plot_data[0].min()

        print "averaging along time axis..."
        plot_max = plot_data.max(axis=0)
        plot_min = plot_data.min(axis=0)
        plot_data = plot_data.mean(axis=0)

        # Rebin to max number of points
        dec_fac_x = 1
        MAX_PLT_POINTS = 8*64  # Low resoluition to see the difference.
        if plot_data.shape[0] > MAX_PLT_POINTS:
            dec_fac_x = plot_data.shape[0] / MAX_PLT_POINTS

        plot_data = rebin(plot_data, dec_fac_x, 1)
        plot_min = rebin(plot_min, dec_fac_x, 1)
        plot_max = rebin(plot_max, dec_fac_x, 1)
        plot_f    = rebin(plot_f, dec_fac_x, 1)

        if logged:
            plt.plot(plot_f, db(plot_data),'k', **kwargs)
            plt.plot(plot_f, db(plot_max),'b', **kwargs)
            plt.plot(plot_f, db(plot_min),'b', **kwargs)
            plt.ylabel("Power [dB]")
        else:
            plt.plot(plot_f, plot_data,'k', **kwargs)
            plt.plot(plot_f, plot_max,'b', **kwargs)
            plt.plot(plot_f, plot_min,'b', **kwargs)
            plt.ylabel("Power [counts]")
        plt.xlabel("Frequency [MHz]")

        try:
            plt.title(self.header['source_name'])
        except KeyError:
            plt.title(self.filename)

        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.xlim(plot_f[0], plot_f[-1])
        plt.ylim(db(fig_min),db(fig_max))

    def plot_waterfall(self, f_start=None, f_stop=None, if_id=0, logged=True,cb=True, **kwargs):
        """ Plot waterfall of data

        Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            logged (bool): Plot in linear (False) or dB units (True),
            cb (bool): for plotting the colorbar
            kwargs: keyword args to be passed to matplotlib imshow()
        """
        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        if logged:
            plot_data = db(plot_data)

        # Make sure waterfall plot is under 4k*4k
        dec_fac_x, dec_fac_y = 1, 1
        if plot_data.shape[0] > MAX_IMSHOW_POINTS[0]:
            dec_fac_x = plot_data.shape[0] / MAX_IMSHOW_POINTS[0]

        if plot_data.shape[1] > MAX_IMSHOW_POINTS[1]:
            dec_fac_y =  plot_data.shape[1] /  MAX_IMSHOW_POINTS[1]

        plot_data = rebin(plot_data, dec_fac_x, dec_fac_y)

        try:
            plt.title(self.header['source_name'])
        except KeyError:
            plt.title(self.filename)

        plt.imshow(plot_data,
            aspect='auto',
            rasterized=True,
            interpolation='nearest',
            extent=(plot_f[0], plot_f[-1], self.timestamps[-1], self.timestamps[0]),
            cmap='viridis',
            **kwargs
        )
        if cb:
            plt.colorbar()
        plt.xlabel("Frequency [MHz]")
        plt.ylabel("Time [MJD]")

    def plot_time_series(self, f_start=None, f_stop=None, if_id=0, logged=True, orientation=None , **kwargs):
        ''' Plot kurtosis

         Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            logged (bool): Plot in linear (False) or dB units (True),
            kwargs: keyword args to be passed to matplotlib imshow()
        '''

        ax = plt.gca()
        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        if logged:
            plot_data = db(plot_data)

        plot_data = plot_data.mean(axis=1)

        if 'v' in orientation:
            plt.plot(plot_data,range(len(plot_data))[::-1], **kwargs)
        else:
            plt.plot(plot_data, **kwargs)
            plt.xlabel("Time [s]")

        ax.autoscale(axis='both',tight=True)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)

    def plot_kurtosis(self, f_start=None, f_stop=None, if_id=0, **kwargs):
        ''' Plot kurtosis

         Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            kwargs: keyword args to be passed to matplotlib imshow()
        '''
        ax = plt.gca()

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)
        plot_kurtossis = np.zeros(len(plot_f))

        for i in range(len(plot_f)):
            plot_kurtossis[i] = scipy.stats.kurtosis(plot_data[:,i],nan_policy='omit')

        plt.plot(plot_f, plot_kurtossis, **kwargs)
        plt.ylabel("Kurtosis")
        plt.xlabel("Frequency [MHz]")

        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.xlim(plot_f[0], plot_f[-1])

    def plot_all(self, t=0, f_start=None, f_stop=None, logged=False, if_id=0, c=None, **kwargs):
        """ Plot waterfall of data as well as spectrum; also, placeholder to make even more complicated plots in the future.

        Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            logged (bool): Plot in linear (False) or dB units (True),
            t (int): integration number to plot (0 -> len(data))
            logged (bool): Plot in linear (False) or dB units (True)
            if_id (int): IF identification (if multiple IF signals in file)
            c: color for line
            kwargs: keyword args to be passed to matplotlib plot() and imshow()
        """

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        nullfmt = NullFormatter()         # no labels

        # definitions for the axes
        left, width = 0.35, 0.5
        bottom, height = 0.45, 0.5
        width2, height2 = 0.1125, 0.15
        bottom2, left2 = bottom-height2-.025, left-width2-.02
        bottom3, left3 = bottom2-height2-.025, 0.075

        rect_waterfall = [left, bottom, width, height]
        rect_colorbar = [left+width, bottom, .025, height]
        rect_spectrum = [left, bottom2, width, height2]
        rect_min_max = [left, bottom3, width, height2]
        rect_timeseries = [left+width, bottom, width2, height]
        rect_kurtosis = [left3, bottom3, 0.25, height2]
        rect_header = [left3-.05, bottom, 0.2, height]

        #--------
        axWaterfall = plt.axes(rect_waterfall)
        print 'Ploting Waterfall'
        self.plot_waterfall(f_start=f_start, f_stop=f_stop,cb=False)
        plt.xlabel('')

        # no labels
        axWaterfall.xaxis.set_major_formatter(nullfmt)

        #--------
#         axColorbar = plt.axes(rect_colorbar)
#         print 'Ploting Colorbar'
#         print plot_data.max()
#         print plot_data.min()
#
#         plot_colorbar = range(plot_data.min(),plot_data.max(),int((plot_data.max()-plot_data.min())/plot_data.shape[0]))
#         plot_colorbar = np.array([[plot_colorbar],[plot_colorbar]])
#
#         plt.imshow(plot_colorbar,aspect='auto', rasterized=True, interpolation='nearest',)

#         axColorbar.xaxis.set_major_formatter(nullfmt)
#         axColorbar.yaxis.set_major_formatter(nullfmt)

#         heatmap = axColorbar.pcolor(plot_data, edgecolors = 'none', picker=True)
#         plt.colorbar(heatmap, cax = axColorbar)

        #--------
        axSpectrum = plt.axes(rect_spectrum)
        print 'Ploting Spectrum'
        self.plot_spectrum(logged=logged, f_start=f_start, f_stop=f_stop, t=t)
        plt.title('')
        axSpectrum.yaxis.tick_right()
        axSpectrum.yaxis.set_label_position("right")
        plt.xlabel('')
        axSpectrum.xaxis.set_major_formatter(nullfmt)

        #--------
        axTimeseries = plt.axes(rect_timeseries)
        print 'Ploting Timeseries'
        self.plot_time_series(f_start=f_start, f_stop=f_stop,orientation='v')
        axTimeseries.yaxis.set_major_formatter(nullfmt)
        axTimeseries.xaxis.set_major_formatter(nullfmt)

        #--------
        axKurtosis = plt.axes(rect_kurtosis)
        print 'Ploting Kurtosis'
        self.plot_kurtosis(f_start=f_start, f_stop=f_stop)

        #--------
        axMinMax = plt.axes(rect_min_max)
        print 'Ploting Min Max'
        self.plot_spectrum_min_max(logged=logged, f_start=f_start, f_stop=f_stop, t=t)
        plt.title('')
        axMinMax.yaxis.tick_right()
        axMinMax.yaxis.set_label_position("right")

        #--------
        axHeader = plt.axes(rect_header)
        print 'Ploting Header'
        plot_header = '\n'.join(['%s:  %s'%(key.upper(),value) for (key,value) in self.header.items() if 'source_name' not in key])
        plt.text(0,1,plot_header,ha='left', va='top', wrap=True)

        axHeader.set_axis_bgcolor('white')
        axHeader.xaxis.set_major_formatter(nullfmt)
        axHeader.yaxis.set_major_formatter(nullfmt)

    def write_to_filterbank(self, filename_out):
        """ Write data to filterbank file.

        Args:
            filename_out (str): Name of output file
        """

        #calibrate data
        #self.data = calibrate(mask(self.data.mean(axis=0)[0]))
        #rewrite header to be consistent with modified data
        self.header['fch1']   = self.freqs[0]
        self.header['foff']   = self.freqs[1] - self.freqs[0]
        self.header['nchans'] = self.freqs.shape[0]
        #self.header['tsamp']  = self.data.shape[0] * self.header['tsamp']

        n_bytes  = self.header['nbits'] / 8
        with open(filename_out, "w") as fileh:
            fileh.write(generate_sigproc_header(self))
            j = self.data
            if n_bytes == 4:
                np.float32(j[:, ::-1].ravel()).tofile(fileh)
            elif n_bytes == 2:
                np.int16(j[:, ::-1].ravel()).tofile(fileh)
            elif n_bytes == 1:
                np.int8(j[:, ::-1].ravel()).tofile(fileh)

    def write_to_hdf5(self, filename_out, *args, **kwargs):
        """ Write data to HDF5 file.

        Args:
            filename_out (str): Name of output file
        """
        if not HAS_HDF5:
            raise RuntimeError("h5py package required for HDF5 output.")

        with h5py.File(filename_out, 'w') as h5:

            dset = h5.create_dataset('data',
                              data=self.data,
                              compression='lzf')

            dset_mask = h5.create_dataset('mask',
                                     shape=self.data.shape,
                                     compression='lzf',
                                     dtype='uint8')

            dset.dims[0].label = "frequency"
            dset.dims[1].label = "feed_id"
            dset.dims[2].label = "time"

            dset_mask.dims[0].label = "frequency"
            dset_mask.dims[1].label = "feed_id"
            dset_mask.dims[2].label = "time"

            # Copy over header information as attributes
            for key, value in self.header.items():
                dset.attrs[key] = value


def cmd_tool(args=None):
    """ Command line tool for plotting and viewing info on filterbank files """    
    
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Command line utility for reading and plotting filterbank files.")

    parser.add_argument('-p', action='store',  default='a', dest='what_to_plot', type=str,
                        help='Show: "w" waterfall (freq vs. time) plot; "s" integrated spectrum plot, \
                             "a" for all available plots and information; and more.')
    parser.add_argument('filename', type=str,
                        help='Name of file to read')
    parser.add_argument('-b', action='store', default=None, dest='f_start', type=float,
                        help='Start frequency (begin), in MHz')
    parser.add_argument('-e', action='store', default=None, dest='f_stop', type=float,
                        help='Stop frequency (end), in MHz')
    parser.add_argument('-B', action='store', default=None, dest='t_start', type=int,
                        help='Start integration (begin) ID')
    parser.add_argument('-E', action='store', default=None, dest='t_stop', type=int,
                        help='Stop integration (end) ID')
    parser.add_argument('-i', action='store_true', default=False, dest='info_only',
                        help='Show info only')
    parser.add_argument('-a', action='store_true', default=False, dest='average',
                       help='average along time axis (plot spectrum only)')
    parser.add_argument('-s', action='store', default='', dest='plt_filename', type=str,
                       help='save plot graphic to file (give filename as argument)')
    parser.add_argument('-S', action='store_true', default=False, dest='save_only',
                       help='Turn off plotting of data and only save to file.')
    parser.add_argument('-D', action='store', default=True, type=int, dest='blank_dc',
                       help='Blank DC bin. Need to know number of coarse channels.')
    args = parser.parse_args()

    # Open filterbank data
    filename = args.filename
    load_data = not args.info_only

    # only load one integration if looking at spectrum
    wtp = args.what_to_plot
    if not wtp or 's' in wtp:
        if args.t_start == None:
            t_start = 0
        else:
            t_start = args.t_start
        t_stop  = t_start + 1

        if args.average:
            t_start = None
            t_stop  = None
    else:
        t_start = args.t_start
        t_stop  = args.t_stop

    fil = Filterbank(filename, f_start=args.f_start, f_stop=args.f_stop,
                     t_start=t_start, t_stop=t_stop,
                     load_data=load_data)
    fil.info()

    # And if we want to plot data, then plot data.
    if not args.info_only:
        # check start & stop frequencies make sense
        #try:
        #    if args.f_start:
        #        print "Start freq: %2.2f" % args.f_start
        #        assert args.f_start >= fil.freqs[0] or np.isclose(args.f_start, fil.freqs[0])
        #
        #    if args.f_stop:
        #        print "Stop freq: %2.2f" % args.f_stop
        #        assert args.f_stop <= fil.freqs[-1] or np.isclose(args.f_stop, fil.freqs[-1])
        #except AssertionError:
        #    print "Error: Start and stop frequencies must lie inside file's frequency range."
        #    print "i.e. between %2.2f-%2.2f MHz." % (fil.freqs[0], fil.freqs[-1])
        #    exit()

        if args.blank_dc:
            print "Blanking DC bin"
            fil.blank_dc(args.blank_dc)

        if "w" in args.what_to_plot:
            plt.figure("waterfall", figsize=(8, 6))
            fil.plot_waterfall(f_start=args.f_start, f_stop=args.f_stop)
        elif "s" in args.what_to_plot:
            plt.figure("Spectrum", figsize=(8, 6))
            fil.plot_spectrum(logged=True, f_start=args.f_start, f_stop=args.f_stop, t='all')
        elif "mm" in args.what_to_plot:
            plt.figure("min max", figsize=(8, 6))
            fil.plot_spectrum_min_max(logged=True, f_start=args.f_start, f_stop=args.f_stop, t='all')
        elif "k" in args.what_to_plot:
            plt.figure("kurtosis", figsize=(8, 6))
            fil.plot_kurtosis(f_start=args.f_start, f_stop=args.f_stop)
        elif "t" in args.what_to_plot:
            plt.figure("Time Series", figsize=(8, 6))
            fil.plot_time_series(f_start=args.f_start, f_stop=args.f_stop)
        elif "a" in args.what_to_plot:
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=args.f_start, f_stop=args.f_stop, t='all')

        if args.plt_filename != '':
            plt.savefig(args.plt_filename)

        if not args.save_only:
            if os.environ.has_key('DISPLAY'):
                plt.show()
            else:
                print "No $DISPLAY available."


if __name__ == "__main__":
    cmd_tool()