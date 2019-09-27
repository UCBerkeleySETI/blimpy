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


import sys

from astropy.time import Time

import logging as logger


from blimpy.utils import unpack_2to8
from blimpy.io.sigproc import *
from blimpy.plotting import *

try:
    import h5py
    HAS_HDF5 = True
except ImportError:
    HAS_HDF5 = False

import warnings

#import pdb #pdb.set_trace()


###
# Config values
###

MAX_DATA_ARRAY_SIZE = 1024 * 1024 * 1024 * 2 # Max size of data array to load into memory
MAX_HEADER_BLOCKS   = 100                    # Max size of header (in 512-byte blocks)


###
# Main blimpy class
###

class Filterbank(object):
    """ Class for loading and plotting blimpy data """

    def __repr__(self):
        return "Filterbank data: %s" % self.filename

    def __init__(self, filename=None, f_start=None, f_stop=None,
                 t_start=None, t_stop=None, load_data=True,
                 header_dict=None, data_array=None, blank_dc=False,cal_band_pass=False):
        """ Class for loading and plotting blimpy data.

        This class parses the blimpy file and stores the header and data
        as objects:
            fb = Filterbank('filename_here.fil')
            fb.header        # blimpy header, as a dictionary
            fb.data          # blimpy data, as a numpy array

        Args:
            filename (str): filename of blimpy file.
            f_start (float): start frequency in MHz
            f_stop (float): stop frequency in MHz
            t_start (int): start integration ID
            t_stop (int): stop integration ID
            load_data (bool): load data. If set to False, only header will be read.
            header_dict (dict): Create blimpy from header dictionary + data array
            data_array (np.array): Create blimpy from header dict + data array
            blank_dc (bool): Blanks the DC bin or not.
        """
        warnings.warn("Filterbank is deprecated. Please use Waterfall()")

        if filename:
            self.filename = filename
            if HAS_HDF5:
                if h5py.is_hdf5(filename):
                    # TODO: self.read_hdf5(filename, f_start, f_stop, t_start, t_stop, load_data)
                    self.read_hdf5(filename, f_start, f_stop, t_start, t_stop, load_data)
                else:
                    self.read_filterbank(filename, f_start, f_stop, t_start, t_stop, load_data)
            else:
                self.read_filterbank(filename, f_start, f_stop, t_start, t_stop, load_data)
        elif header_dict is not None and data_array is not None:
            self.filename = b''
            self.header = header_dict
            self.data = data_array
            self.n_ints_in_file = 0

            self._setup_freqs()

        else:
            pass

    def read_hdf5(self, filename, f_start=None, f_stop=None,
                        t_start=None, t_stop=None, load_data=True):
        """ Populate Filterbank instance with data from HDF5 file

        Note:
            This is to be deprecated in future, please use Waterfall() to open files.
        """
        print("Warning: this function will be deprecated in the future. Please use Waterfall to open HDF5 files.")
#        raise DeprecationWarning('Please use Waterfall to open HDF5 files.')

        self.header = {}
        self.filename = filename
        self.h5 = h5py.File(filename)
        for key, val in self.h5[b'data'].attrs.items():
            if six.PY3:
                key = bytes(key, 'ascii')
            if key == b'src_raj':
                self.header[key] = Angle(val, unit='hr')
            elif key == b'src_dej':
                self.header[key] = Angle(val, unit='deg')
            else:
                self.header[key] = val

        self.n_ints_in_file = self.h5[b"data"].shape[0]
        i_start, i_stop, chan_start_idx, chan_stop_idx = self._setup_freqs(f_start=f_start, f_stop=f_stop)
        ii_start, ii_stop, n_ints = self._setup_time_axis(t_start=t_start, t_stop=t_stop)

        if load_data:
            self.data = self.h5[b"data"][ii_start:ii_stop, :, chan_start_idx:chan_stop_idx]

            self.file_size_bytes = os.path.getsize(self.filename)

#         if self.header[b'foff'] < 0:
#             self.data = self.data[..., ::-1] # Reverse data

        else:
            print("Skipping data load...")
            self.data = np.array([0])
            self.n_ints_in_file  = 0
            self.file_size_bytes = os.path.getsize(self.filename)


    def _setup_freqs(self, f_start=None, f_stop=None):
        """ Setup frequency axis """
        ## Setup frequency axis
        f0 = self.header[b'fch1']
        f_delt = self.header[b'foff']

        i_start, i_stop = 0, self.header[b'nchans']
        if f_start:
            i_start = int((f_start - f0) / f_delt)
        if f_stop:
            i_stop  = int((f_stop - f0)  / f_delt)

        #calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx  = np.int(i_stop)

        #create freq array
        if i_start < i_stop:
            i_vals = np.arange(chan_start_idx, chan_stop_idx)
        else:
            i_vals = np.arange(chan_stop_idx, chan_start_idx)

        self.freqs = f_delt * i_vals + f0

#         if f_delt < 0:
#             self.freqs = self.freqs[::-1]

        if chan_stop_idx < chan_start_idx:
            chan_stop_idx, chan_start_idx = chan_start_idx,chan_stop_idx

        return i_start, i_stop, chan_start_idx, chan_stop_idx

    def _setup_time_axis(self, t_start=None, t_stop=None):
        """  Setup time axis. """

        # now check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if t_start:
            ii_start = t_start
        if t_stop:
            ii_stop = t_stop
        n_ints = ii_stop - ii_start

        ## Setup time axis
        t0 = self.header[b'tstart']
        t_delt = self.header[b'tsamp']

        self.timestamps = np.arange(0, n_ints) * t_delt / 24./60./60 + t0

        return ii_start, ii_stop, n_ints

    def read_filterbank(self, filename=None, f_start=None, f_stop=None,
                        t_start=None, t_stop=None, load_data=True):
        """ Populate Filterbank instance with data from Filterbank file

        Note:
            This is to be deprecated in future, please use Waterfall() to open files.
        """
        if filename is None:
            filename = self.filename
        else:
            self.filename = filename

        self.header = read_header(filename)

        #convert input frequencies into what their corresponding index would be
        i_start, i_stop, chan_start_idx, chan_stop_idx = self._setup_freqs(f_start=f_start, f_stop=f_stop)

        n_bits  = self.header[b'nbits']
        n_bytes  = int(self.header[b'nbits'] / 8)
        n_chans = self.header[b'nchans']
        n_chans_selected = self.freqs.shape[0]
        n_ifs   = self.header[b'nifs']

        # Load binary data
        self.idx_data = len_header(filename)
        f = open(filename, 'rb')
        f.seek(self.idx_data)
        filesize = os.path.getsize(self.filename)
        n_bytes_data = filesize - self.idx_data

        # Finally add some other info to the class as objects
        self.n_ints_in_file  = calc_n_ints_in_file(self.filename)
        self.file_size_bytes = filesize

        ## Setup time axis
        ii_start, ii_stop, n_ints = self._setup_time_axis(t_start=t_start, t_stop=t_stop)

        # Seek to first integration
        f.seek(int(ii_start * n_bits * n_ifs * n_chans / 8), 1)

        # Set up indexes used in file read (taken out of loop for speed)
        i0 = np.min((chan_start_idx, chan_stop_idx))
        i1 = np.max((chan_start_idx, chan_stop_idx))

        #Set up the data type (taken out of loop for speed)
        if n_bits == 2:
            dd_type = b'uint8'
            n_chans_selected = int(n_chans_selected/4)
        elif n_bytes == 4:
            dd_type = b'float32'
        elif n_bytes == 2:
            dd_type = b'uint16'
        elif n_bytes == 1:
            dd_type = b'uint8'

        if load_data:

            if n_ints * n_ifs * n_chans_selected > MAX_DATA_ARRAY_SIZE:
                print("[Filterbank]  Error: data array is too large to load. Either select fewer points or manually increase MAX_DATA_ARRAY_SIZE. Large files are now handle with Waterfall .")
                sys.exit()

            if n_bits == 2:
                self.data = np.zeros((n_ints, n_ifs, n_chans_selected*4), dtype=dd_type)
            else:
                self.data = np.zeros((n_ints, n_ifs, n_chans_selected), dtype=dd_type)

            for ii in range(n_ints):
                """d = f.read(n_bytes * n_chans * n_ifs)
                """

                for jj in range(n_ifs):

                    f.seek(n_bytes * i0, 1) # 1 = from current location
                    #d = f.read(n_bytes * n_chans_selected)
                    #bytes_to_read = n_bytes * n_chans_selected

                    dd = np.fromfile(f, count=n_chans_selected, dtype=dd_type)

                    # Reverse array if frequency axis is flipped
#                     if f_delt < 0:
#                         dd = dd[::-1]

                    if n_bits == 2:
                        dd = unpack_2to8(dd)
                    self.data[ii, jj] = dd

                    f.seek(n_bytes * (n_chans - i1), 1)  # Seek to start of next block
        else:
            print("Skipping data load...")
            self.data = np.array([0], dtype=dd_type)

    def info(self):
        """ Print header information """

        for key, val in self.header.items():
            if key == b'src_raj':
                val = val.to_string(unit=u.hour, sep=':')
            if key == b'src_dej':
                val = val.to_string(unit=u.deg, sep=':')
            if key == b'tsamp':
                val *= u.second
            if key in ('foff', 'fch1'):
                val *= u.MHz
            if key == b'tstart':
                print("%16s : %32s" % ("tstart (ISOT)", Time(val, format='mjd').isot))
                key = "tstart (MJD)"
            print("%16s : %32s" % (key, val))

        print("\n%16s : %32s" % ("Num ints in file", self.n_ints_in_file))
        print("%16s : %32s" % ("Data shape", self.data.shape))
        print("%16s : %32s" % ("Start freq (MHz)", self.freqs[0]))
        print("%16s : %32s" % ("Stop freq (MHz)", self.freqs[-1]))

    def grab_data(self, f_start=None, f_stop=None, t_start=None, t_stop=None, if_id=0):
        """ Extract a portion of data by frequency range.

        Args:
            f_start (float): start frequency in MHz
            f_stop (float): stop frequency in MHz
            if_id (int): IF input identification (req. when multiple IFs in file)

        Returns:
            (freqs, data) (np.arrays): frequency axis in MHz and data subset
        """

        if f_start is None:
            f_start = self.freqs[0]
        if f_stop is None:
            f_stop = self.freqs[-1]

        i0 = np.argmin(np.abs(self.freqs - f_start))
        i1 = np.argmin(np.abs(self.freqs - f_stop))

        if i0 < i1:
            plot_f    = self.freqs[i0:i1 + 1]
            plot_data = np.squeeze(self.data[t_start:t_stop, ..., i0:i1 + 1])
        else:
            plot_f    = self.freqs[i1:i0 + 1]
            plot_data = np.squeeze(self.data[t_start:t_stop, ..., i1:i0 + 1])

        return plot_f, plot_data

    def write_to_fil(self, filename_out):
        """ Write data to blimpy file.

        Args:
            filename_out (str): Name of output file
        """

        print("[Filterbank] Warning: Non-standard function to write in filterbank (.fil) format. Please use Waterfall.")

        n_bytes  = int(self.header[b'nbits'] / 8)
        with open(filename_out, "wb") as fileh:
            fileh.write(generate_sigproc_header(self))
            j = self.data
            if n_bytes == 4:
                np.float32(j.ravel()).tofile(fileh)
            elif n_bytes == 2:
                np.int16(j.ravel()).tofile(fileh)
            elif n_bytes == 1:
                np.int8(j.ravel()).tofile(fileh)

    def write_to_hdf5(self, filename_out, *args, **kwargs):
        """ Write data to HDF5 file.

        Args:
            filename_out (str): Name of output file
        """

        print("[Filterbank] Warning: Non-standard function to write in HDF5 (.h5) format. Please use Waterfall.")

        if not HAS_HDF5:
            raise RuntimeError("h5py package required for HDF5 output.")

        with h5py.File(filename_out, 'w') as h5:

            dset = h5.create_dataset(b'data',
                              data=self.data,
                              compression='lzf')

            dset_mask = h5.create_dataset(b'mask',
                                     shape=self.data.shape,
                                     compression='lzf',
                                     dtype='uint8')

            dset.dims[0].label = b"frequency"
            dset.dims[1].label = b"feed_id"
            dset.dims[2].label = b"time"

            dset_mask.dims[0].label = b"frequency"
            dset_mask.dims[1].label = b"feed_id"
            dset_mask.dims[2].label = b"time"

            # Copy over header information as attributes
            for key, value in self.header.items():
                dset.attrs[key] = value

    def plot_spectrum(self, t=0, f_start=None, f_stop=None, logged=False, if_id=0, c=None, **kwargs):
        plot_spectrum(self, t, f_start, f_stop, logged, if_id, c, **kwargs)

    def plot_waterfall(self, f_start=None, f_stop=None, if_id=0, logged=True, cb=True, MJD_time=False, **kwargs):
        plot_waterfall(self, f_start, f_stop, if_id, logged, cb, MJD_time, **kwargs)

    def plot_time_series(self, f_start=None, f_stop=None, if_id=0, logged=True, orientation='h', MJD_time=False,
                         **kwargs):
        plot_time_series(self, f_start, f_stop, if_id, logged, orientation, MJD_time, **kwargs)

    def plot_spectrum_min_max(self, t=0, f_start=None, f_stop=None, logged=False, if_id=0, c=None, **kwargs):
        plot_spectrum_min_max(self, t, f_start, f_stop, logged, if_id, c, **kwargs)

    def plot_kurtosis(self, f_start=None, f_stop=None, if_id=0, **kwargs):
        plot_kurtosis(self, f_start, f_stop, if_id, **kwargs)

    def plot_all(self, t=0, f_start=None, f_stop=None, logged=False, if_id=0, kurtosis=True, **kwargs):
        plot_all(self, t, f_start, f_stop, logged, if_id, kurtosis, **kwargs)


def cmd_tool(args=None):
    """ Command line tool for plotting and viewing info on filterbank files """

    from argparse import ArgumentParser
    from plotting.config import plt

    parser = ArgumentParser(description="Command line utility for reading and plotting filterbank files.")

    parser.add_argument('-p', action='store',  default='ank', dest='what_to_plot', type=str,
                        help='Show: "w" waterfall (freq vs. time) plot; "s" integrated spectrum plot; \
                        "t" for time series; "mm" for spectrum including min max; "k" for kurtosis; \
                        "a" for all available plots and information; and "ank" for all but kurtosis.')
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
    parser.add_argument('-D', action='store_false', default=True, dest='blank_dc',
                       help='Use to not blank DC bin.')
    parser.add_argument('-c', action='store_true', default=False, dest='calibrate_band_pass',
                        help='Calibrate band pass.')
    args = parser.parse_args()

    # Open blimpy data
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

    if args.info_only:
        args.blank_dc = False
        args.calibrate_band_pass = False

    fil = Filterbank(filename, f_start=args.f_start, f_stop=args.f_stop,
                     t_start=t_start, t_stop=t_stop,
                     load_data=load_data,blank_dc=args.blank_dc,
                     cal_band_pass=args.calibrate_band_pass)
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

        if args.what_to_plot == "w":
            plt.figure("waterfall", figsize=(8, 6))
            fil.plot_waterfall(f_start=args.f_start, f_stop=args.f_stop)
        elif args.what_to_plot == "s":
            plt.figure("Spectrum", figsize=(8, 6))
            fil.plot_spectrum(logged=True, f_start=args.f_start, f_stop=args.f_stop, t='all')
        elif args.what_to_plot == "mm":
            plt.figure("min max", figsize=(8, 6))
            fil.plot_spectrum_min_max(logged=True, f_start=args.f_start, f_stop=args.f_stop, t='all')
        elif args.what_to_plot == "k":
            plt.figure("kurtosis", figsize=(8, 6))
            fil.plot_kurtosis(f_start=args.f_start, f_stop=args.f_stop)
        elif args.what_to_plot == "t":
            plt.figure("Time Series", figsize=(8, 6))
            fil.plot_time_series(f_start=args.f_start, f_stop=args.f_stop,orientation='h')
        elif args.what_to_plot == "a":
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=args.f_start, f_stop=args.f_stop, t='all')
        elif args.what_to_plot == "ank":
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=args.f_start, f_stop=args.f_stop, t='all',kurtosis=False)

        if args.plt_filename != '':
            plt.savefig(args.plt_filename)

        if not args.save_only:
            if 'DISPLAY' in os.environ.keys():
                plt.show()
            else:
                print("No $DISPLAY available.")


if __name__ == "__main__":
    cmd_tool()
