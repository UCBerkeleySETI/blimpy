#!/usr/bin/env python
"""
# filterbank.py

Python class and command line utility for reading and plotting filterbank files.

This provides a class, Filterbank(), which can be used to read a .fil file:


    fil = Filterbank('test_psr.fil')
    print(fil.header)
    print(fil.data.shape)
    print(fil.freqs)

    plt.figure()
    fil.plot_spectrum(t=0)
    plt.show()


TODO: check the file seek logic works correctly for multiple IFs

"""


import sys
import six

from astropy.time import Time
import scipy.stats
import numpy as np
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter

import logging as logger

try:
    from .utils import db, lin, rebin, closest, unpack_2to8
    from .sigproc import *
except:
    from utils import db, lin, rebin, closest, unpack_2to8
    from sigproc import *

try:
    import h5py
    HAS_HDF5 = True
except ImportError:
    HAS_HDF5 = False

try:
    from pyslalib import slalib as s
    HAS_SLALIB = True
except ImportError:
    HAS_SLALIB = False

#import pdb #pdb.set_trace()

# Check if $DISPLAY is set (for handling plotting on remote machines with no X-forwarding)
if 'DISPLAY' in os.environ.keys():
    import pylab as plt
else:
    import matplotlib
    matplotlib.use('Agg')
    import pylab as plt

from matplotlib import gridspec
plt.rcParams['axes.formatter.useoffset'] = False


###
# Config values
###

MAX_PLT_POINTS      = 65536                  # Max number of points in matplotlib plot
MAX_IMSHOW_POINTS   = (8192, 4096)           # Max number of points in imshow plot
MAX_DATA_ARRAY_SIZE = 1024 * 1024 * 1024 * 2 # Max size of data array to load into memory
MAX_HEADER_BLOCKS   = 100                    # Max size of header (in 512-byte blocks)


# Telescope coordinates (needed for LSR calc)
parkes_coords = (-32.998370, 148.263659,  324.00)
gbt_coords    = (38.4331294, 79.8398397, 824.36)


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

        if blank_dc:
            print("Blanking DC bin")
            n_coarse_chan = self.calc_n_coarse_chan()
            self.blank_dc(n_coarse_chan)

        if cal_band_pass:
            print("Calibrating the band pass.")
            self.calibrate_band_pass_N1()

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

    def compute_lst(self):
        """ Compute LST for observation """
        if self.header[b'telescope_id'] == 6:
            self.coords = gbt_coords
        elif self.header[b'telescope_id'] == 4:
            self.coords = parkes_coords
        else:
            raise RuntimeError("Currently only Parkes and GBT supported")
        if HAS_SLALIB:
            # dut1 = (0.2 /3600.0) * np.pi/12.0
            dut1 = 0.0
            mjd = self.header[b'tstart']
            tellong = np.deg2rad(self.coords[1])
            last = s.sla_gmst(mjd) - tellong + s.sla_eqeqx(mjd) + dut1
            # lmst = s.sla_gmst(mjd) - tellong
            if last < 0.0 : last = last + 2.0*np.pi
            return last
        else:
            raise RuntimeError("This method requires pySLALIB")

    def compute_lsrk(self):
        """ Computes the LSR in km/s

        uses the MJD, RA and DEC of observation to compute
        along with the telescope location. Requires pyslalib
        """
        ra = Angle(self.header[b'src_raj'], unit='hourangle')
        dec = Angle(self.header[b'src_dej'], unit='degree')
        mjdd = self.header[b'tstart']
        rarad = ra.to('radian').value
        dcrad = dec.to('radian').value
        last = self.compute_lst()
        tellat  = np.deg2rad(self.coords[0])
        tellong = np.deg2rad(self.coords[1])

        # convert star position to vector
        starvect = s.sla_dcs2c(rarad, dcrad)

        # velocity component in ra,dec due to Earth rotation
        Rgeo = s.sla_rverot( tellat, rarad, dcrad, last)

        # get Barycentric and heliocentric velocity and position of the Earth.
        evp = s.sla_evp(mjdd, 2000.0)
        dvb = evp[0]   # barycentric velocity vector, in AU/sec
        dpb = evp[1]   # barycentric position vector, in AU
        dvh = evp[2]   # heliocentric velocity vector, in AU/sec
        dph = evp[3]   # heliocentric position vector, in AU

        # dot product of vector to object and heliocentric velocity
        # convert AU/sec to km/sec
        vcorhelio = -s.sla_dvdv( starvect, dvh) *149.597870e6
        vcorbary  = -s.sla_dvdv( starvect, dvb) *149.597870e6

        # rvlsrd is velocity component in ra,dec direction due to the Sun's
        # motion with respect to the "dynamical" local standard of rest
        rvlsrd = s.sla_rvlsrd( rarad,dcrad)

        # rvlsrk is velocity component in ra,dec direction due to i
        # the Sun's motion w.r.t the "kinematic" local standard of rest
        rvlsrk = s.sla_rvlsrk( rarad,dcrad)

        # rvgalc is velocity component in ra,dec direction due to
        #the rotation of the Galaxy.
        rvgalc = s.sla_rvgalc( rarad,dcrad)
        totalhelio = Rgeo + vcorhelio
        totalbary  = Rgeo + vcorbary
        totallsrk = totalhelio + rvlsrk
        totalgal  = totalbary  + rvlsrd + rvgalc

        return totallsrk

    def blank_dc(self, n_coarse_chan):
        """ Blank DC bins in coarse channels.

        Note: currently only works if entire file is read, and in GBT data.

        """

        if self.header[b'telescope_id'] != 6:
            logger.info('Data in not from GBT, ignoring blanking of DC channel.')
            return None

        if n_coarse_chan < 1:
            logger.warning('Coarse channel number < 1, unable to blank DC channel.')
            return None

        if not n_coarse_chan % int(n_coarse_chan) == 0:
            logger.warning('Selection does not contain an integer number of coarse channels, unable to blank DC channel.')
            return None

        # n_coarse_chan is already integer, but now using it with integer format.
        n_coarse_chan = int(n_coarse_chan)

        n_chan = self.data.shape[-1]
        n_chan_per_coarse = n_chan // n_coarse_chan
        mid_chan = n_chan_per_coarse // 2

        for ii in range(n_coarse_chan):
            ss = ii*n_chan_per_coarse
            self.data[..., ss+mid_chan] = np.median(self.data[..., ss+mid_chan+5:ss+mid_chan+10])

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

    def generate_freqs(self, f_start, f_stop):
        """
        returns frequency array [f_start...f_stop]
        """

        fch1 = self.header[b'fch1']
        foff = self.header[b'foff']

        #convert input frequencies into what their corresponding index would be
        i_start = int((f_start - fch1) / foff)
        i_stop  = int((f_stop - fch1)  / foff)

        #calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx  = np.int(i_stop)

        #create freq array
        i_vals = np.arange(chan_stop_idx, chan_start_idx, 1)

        freqs = foff * i_vals + fch1

        return freqs

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

        if plot_data.dtype == 'uint8':
            plot_data = plot_data.astype('uint16')

        return plot_f, plot_data

    def _calc_extent(self,plot_f=None,plot_t=None,MJD_time=False,inv_t=False):
        """ Setup waterfall plotting edges.
        """

        plot_f_begin = plot_f[0]
        plot_f_end = plot_f[-1] + (plot_f[1]-plot_f[0])

        plot_t_begin = self.timestamps[0]
        plot_t_end  = self.timestamps[-1] + (self.timestamps[1] - self.timestamps[0])

        #Returns time in seconds from start of scan.
        if MJD_time == False:
            plot_t_end = (plot_t_end-plot_t_begin)*24.*60.*60
            plot_t_begin = 0.0

        if inv_t == True:
            extent=(plot_f_begin, plot_f_end, plot_t_end, plot_t_begin)
        else:
            extent=(plot_f_begin, plot_f_end, plot_t_begin, plot_t_end)

        return extent

    def plot_spectrum(self, t=0, f_start=None, f_stop=None, logged=False, tight=False, if_id=0, c=None, **kwargs):
        """ Plot frequency spectrum of a given file

        Args:
            t (int): integration number to plot (0 -> len(data))
            logged (bool): Plot in linear (False) or dB units (True)
            tight (bool): If True, restrict intensity labels to fit within bounds of plot. Otherwise, use matplotlib default labels.
            if_id (int): IF identification (if multiple IF signals in file)
            c: color for line
            kwargs: keyword args to be passed to matplotlib plot()
        """
        if self.header[b'nbits'] <=2:
            logged = False
            t='all'
        ax = plt.gca()

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        #Using accending frequency for all plots.
        if self.header[b'foff'] < 0:
            plot_data = plot_data[..., ::-1] # Reverse data
            plot_f = plot_f[::-1]

        if isinstance(t, int):
            print("extracting integration %i..." % t)
            plot_data = plot_data[t]
        elif t == b'all':
            print("averaging along time axis...")
            #Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
            if len(plot_data.shape) > 1:
                plot_data = plot_data.mean(axis=0)
            else:
                plot_data = plot_data.mean()
        else:
            raise RuntimeError("Unknown integration %s" % t)

        # Rebin to max number of points
        dec_fac_x = 1
        if plot_data.shape[0] > MAX_PLT_POINTS:
            dec_fac_x = int(plot_data.shape[0] / MAX_PLT_POINTS)

        plot_data = rebin(plot_data, dec_fac_x, 1)
        plot_f    = rebin(plot_f, dec_fac_x, 1)

        if not c:
            kwargs['c'] = '#333333'


        db_y_min = db(min(plot_data))
        db_y_max = db(max(plot_data))
        if not logged and not tight:
            pass
        else:
            if logged:
                y_min = db_y_min
                y_max = db_y_max
            else:
                y_min = min(plot_data)
                y_max = max(plot_data)

                order_of_mag = int(np.floor(np.log10(y_max)))

                y_min /= 10**order_of_mag
                y_max /= 10**order_of_mag
            spread = y_max - y_min
            plt.ylim(y_min - spread / 4, y_max + spread / 4)

        if logged:
            plt.plot(plot_f, db(plot_data),label='Stokes I', **kwargs)
            plt.ylabel("Power [dB]")
        elif tight:
            plt.plot(plot_f, plot_data / 10**order_of_mag,label='Stokes I', **kwargs)
            plt.ylabel(r"Power [cts / $10^{%d}$]" % (order_of_mag))
        else:
            plt.plot(plot_f, plot_data,label='Stokes I', **kwargs)
            plt.ylabel("Power [counts]")
        plt.xlabel("Frequency [MHz]")
        plt.legend()

        try:
            plt.title(self.header[b'source_name'])
        except KeyError:
            plt.title(self.filename)

        plt.xlim(plot_f[0], plot_f[-1])

    def plot_spectrum_min_max(self, t=0, f_start=None, f_stop=None, logged=False, tight=False, if_id=0, c=None, **kwargs):
        """ Plot frequency spectrum of a given file

        Args:
            logged (bool): Plot in linear (False) or dB units (True)
            tight (bool): If True, restrict intensity labels to fit within bounds of plot. Otherwise, use matplotlib default labels.
            if_id (int): IF identification (if multiple IF signals in file)
            c: color for line
            kwargs: keyword args to be passed to matplotlib plot()
        """
        ax = plt.gca()

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        #Using accending frequency for all plots.
        if self.header[b'foff'] < 0:
            plot_data = plot_data[..., ::-1] # Reverse data
            plot_f = plot_f[::-1]

        ########
        # Why do this?
        fig_max = plot_data[0].max()
        fig_min = plot_data[0].min()
        ########

        print("averaging along time axis...")

        #Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
        if len(plot_data.shape) > 1:
            plot_max = plot_data.max(axis=0)
            plot_min = plot_data.min(axis=0)
            plot_data = plot_data.mean(axis=0)
        else:
            plot_max = plot_data.max()
            plot_min = plot_data.min()
            plot_data = plot_data.mean()

        # Rebin to max number of points
        dec_fac_x = 1
        MAX_PLT_POINTS = 8*64  # Low resoluition to see the difference.
        if plot_data.shape[0] > MAX_PLT_POINTS:
            dec_fac_x = int(plot_data.shape[0] / MAX_PLT_POINTS)

        plot_data = rebin(plot_data, dec_fac_x, 1)
        plot_min = rebin(plot_min, dec_fac_x, 1)
        plot_max = rebin(plot_max, dec_fac_x, 1)
        plot_f    = rebin(plot_f, dec_fac_x, 1)



        db_y_min = db(min(plot_data))
        db_y_max = db(max(plot_data))
        if not logged and not tight:
            pass
        else:
            if logged:
                y_min = db_y_min
                y_max = db_y_max
            elif tight:
                y_min = min(plot_data)
                y_max = max(plot_data)

                order_of_mag = int(np.floor(np.log10(y_max)))

                y_min /= 10**order_of_mag
                y_max /= 10**order_of_mag
            spread = y_max - y_min
            plt.ylim(y_min - spread / 4, y_max + spread / 4)

        if logged:
            plt.plot(plot_f, db(plot_data), "#333333", label='mean', **kwargs)
            plt.plot(plot_f, db(plot_max),  "#e74c3c", label='max', **kwargs)
            plt.plot(plot_f, db(plot_min),  '#3b5b92', label='min', **kwargs)
            plt.ylabel("Power [dB]")
        elif tight:
            plt.plot(plot_f, plot_data / 10**order_of_mag,  "#333333", label='mean', **kwargs)
            plt.plot(plot_f, plot_max / 10**order_of_mag,   "#e74c3c", label='max', **kwargs)
            plt.plot(plot_f, plot_min / 10**order_of_mag,   '#3b5b92', label='min', **kwargs)
            plt.ylabel(r"Power [cts / $10^{%d}$]" % (order_of_mag))
        else:
            plt.plot(plot_f, plot_data,  "#333333", label='mean', **kwargs)
            plt.plot(plot_f, plot_max,   "#e74c3c", label='max', **kwargs)
            plt.plot(plot_f, plot_min,   '#3b5b92', label='min', **kwargs)
            plt.ylabel("Power [counts]")
        plt.xlabel("Frequency [MHz]")
        plt.legend()

        try:
            plt.title(self.header[b'source_name'])
        except KeyError:
            plt.title(self.filename)

        plt.xlim(plot_f[0], plot_f[-1])

    def plot_waterfall(self, f_start=None, f_stop=None, if_id=0, logged=True, cb=True, MJD_time=False, **kwargs):
        """ Plot waterfall of data

        Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            logged (bool): Plot in linear (False) or dB units (True),
            cb (bool): for plotting the colorbar
            kwargs: keyword args to be passed to matplotlib imshow()
        """

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        #Using accending frequency for all plots.
        if self.header[b'foff'] < 0:
            plot_data = plot_data[..., ::-1] # Reverse data
            plot_f = plot_f[::-1]

        if logged == True:
            plot_data = db(plot_data)

        # Make sure waterfall plot is under 4k*4k
        dec_fac_x, dec_fac_y = 1, 1
        if plot_data.shape[0] > MAX_IMSHOW_POINTS[0]:
            dec_fac_x = int(plot_data.shape[0] / MAX_IMSHOW_POINTS[0])

        if plot_data.shape[1] > MAX_IMSHOW_POINTS[1]:
            dec_fac_y =  int(plot_data.shape[1] /  MAX_IMSHOW_POINTS[1])

        plot_data = rebin(plot_data, dec_fac_x, dec_fac_y)

        try:
            plt.title(self.header[b'source_name'])
        except KeyError:
            plt.title(self.filename)

        plt_args = {
            'aspect':'auto',
            'origin':'lower',
            'rasterized':True,
            'interpolation':'nearest',
            'cmap':'viridis'
            }

        for k, v in kwargs.items():
            plt_args[k] = v

        if plt_args['origin'] == 'upper':
            extent = self._calc_extent(plot_f=plot_f,plot_t=self.timestamps,MJD_time=MJD_time,inv_t=True)
        elif plt_args['origin'] == 'lower':
            extent = self._calc_extent(plot_f=plot_f,plot_t=self.timestamps,MJD_time=MJD_time)
        else:
            raise ValueError('The origin parameter not well set:',plt_args['origin'] )

        plt_args['extent'] = extent

        plt.imshow(plot_data, **plt_args)

        if cb:
            plt.colorbar()
        plt.xlabel("Frequency [MHz]")
        if MJD_time:
            plt.ylabel("Time [MJD]")
        else:
            plt.ylabel("Time [s]")

    def plot_time_series(self, f_start=None, f_stop=None, if_id=0, logged=True, tight=False, orientation='h', MJD_time=False, **kwargs):
        """ Plot the time series.

         Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            logged (bool): Plot in linear (False) or dB units (True)
            tight (bool): If True, restrict intensity labels to fit within bounds of plot. Otherwise, use matplotlib default labels.
            kwargs: keyword args to be passed to matplotlib imshow()
        """

        ax = plt.gca()
        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)


        #Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
        if len(plot_data.shape) > 1:
            plot_data = plot_data.mean(axis=1)
        else:
            plot_data = plot_data.mean()

        db_y_min = db(min(plot_data))
        db_y_max = db(max(plot_data))
        if not logged and not tight:
            pass
        else:
            if logged:
                y_min = db_y_min
                y_max = db_y_max
            else:
                y_min = min(plot_data)
                y_max = max(plot_data)

                order_of_mag = int(np.floor(np.log10(y_max)))

                y_min /= 10**order_of_mag
                y_max /= 10**order_of_mag

                plot_data /= 10**order_of_mag 
            spread = y_max - y_min

        if logged and self.header[b'nbits'] >= 8:
            plot_data = db(plot_data)


        #Make proper time axis for plotting (but only for plotting!). Note that this makes the values inclusive.
        extent = self._calc_extent(plot_f=plot_f,plot_t=self.timestamps,MJD_time=MJD_time)
        plot_t = np.linspace(extent[2],extent[3],len(self.timestamps))

        if MJD_time:
            tlabel = "Time [MJD]"
        else:
            tlabel = "Time [s]"

        if logged:
            plabel = "Power [dB]"
        elif tight: 
            plabel = r"Power [cts / $10^{%d}$]" % (order_of_mag)
        else:
            plabel = "Power [counts]"

        # Reverse oder if vertical orientation.
        if 'v' in orientation:
            plt.plot(plot_data, plot_t, **kwargs)
            plt.xlabel(plabel)
            if tight:
                plt.xlim(y_min - spread / 4, y_max + spread / 4)
                plt.xticks([y_min + spread / 4, y_max - spread / 4, y_max + spread / 4])
                ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

                plt.ylim(plot_t[0], plot_t[-1])

        else:
            plt.plot(plot_t, plot_data, **kwargs)
            plt.xlabel(tlabel)
            plt.ylabel(plabel)
            if tight:
                plt.xlim(plot_t[0], plot_t[-1])
                plt.yticks([y_min + spread / 4, y_max - spread / 4, y_max + spread / 4])
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
        if not tight:
            ax.autoscale(axis='both',tight=True)
            
    def plot_kurtosis(self, f_start=None, f_stop=None, if_id=0, **kwargs):
        """ Plot kurtosis
         Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            kwargs: keyword args to be passed to matplotlib imshow()
        """
        ax = plt.gca()

        plot_f, plot_data = self.grab_data(f_start, f_stop, if_id)

        #Using accending frequency for all plots.
        if self.header[b'foff'] < 0:
            plot_data = plot_data[..., ::-1] # Reverse data
            plot_f = plot_f[::-1]

        try:
            plot_kurtosis = scipy.stats.kurtosis(plot_data, axis=0, nan_policy='omit')
        except:
            plot_kurtosis = plot_data*0.0

        plt.plot(plot_f, plot_kurtosis, **kwargs)
        plt.ylabel("Kurtosis")
        plt.xlabel("Frequency [MHz]")

        plt.xlim(plot_f[0], plot_f[-1])

    def plot_all(self, t=0, f_start=None, f_stop=None, logged=False, if_id=0, kurtosis=True, custom_spacing=False, hspace=0., wspace=0., **kwargs):
        """ Plot waterfall of data as well as spectrum; also, placeholder to make even more complicated plots in the future.

        Args:
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            logged (bool): Plot in linear (False) or dB units (True),
            t (int): integration number to plot (0 -> len(data))
            logged (bool): Plot in linear (False) or dB units (True)
            if_id (int): IF identification (if multiple IF signals in file)
            custom_spacing (bool): If True, allows customization of horizontal and vertical spacing via wspace and hspace.
            hspace (float): Controls vertical (height) spacing between plot components in plt.GridSpec.
            wspace (float): Controls horizontal (width) spacing between plot components in plt.GridSpec.
            kwargs: keyword args to be passed to matplotlib plot() and imshow()
        """
        try:
            nbits = self.header[b'nbits']
        except KeyError:
            nbits = self.header['nbits']
        if nbits <=2:
            logged = False

        nullfmt = NullFormatter()  # no labels

        gs = gridspec.GridSpec(24, 32)

        if custom_spacing:
            gs.update(hspace=hspace, wspace=wspace)
        else:
            gs.update(hspace=0, wspace=0)

        # --------
        axMinMax = plt.subplot(gs[20:24, 12:-4])
    #         axMinMax = plt.axes(rect_min_max)
        print('Plotting Min Max')
        self.plot_spectrum_min_max(logged=logged, tight=True, f_start=f_start, f_stop=f_stop, t=t)
        plt.title('')
        axMinMax.yaxis.tick_left()
        axMinMax.yaxis.set_label_position("left")


        axMinMax.xaxis.grid(True, linestyle='--')

        # --------
        axSpectrum = plt.subplot(gs[16:20, 12:-4], sharex=axMinMax)
    #         axSpectrum = plt.axes(rect_spectrum,sharex=axMinMax)
        print('Plotting Spectrum')
        self.plot_spectrum(logged=logged, tight=True, f_start=f_start, f_stop=f_stop, t=t)
        plt.title('')
        axSpectrum.yaxis.tick_left()
        axSpectrum.yaxis.set_label_position("left")
        axSpectrum.set_xlabel('')
    #        axSpectrum.xaxis.set_major_formatter(nullfmt)
        plt.setp(axSpectrum.get_xticklabels(), visible=False)

        axSpectrum.xaxis.grid(True, linestyle='--')

        # --------
        axWaterfall = plt.subplot(gs[4:16, 12:-4], sharex=axMinMax)
    #         axWaterfall = plt.axes(rect_waterfall,sharex=axMinMax)
        print('Plotting Waterfall')
        self.plot_waterfall(f_start=f_start, f_stop=f_stop, logged=logged, cb=False)
        plt.xlabel('')

        # no labels
    #        axWaterfall.xaxis.set_major_formatter(nullfmt)
        plt.setp(axWaterfall.get_xticklabels(), visible=False)

        cax = plt.subplot(gs[3:4, 12:-4])
        plt.colorbar(orientation='horizontal', cax=cax)
        cax.xaxis.set_ticks_position('top') 
        cax.xaxis.set_label_position('top') 
        if logged:
            cax.set_xlabel('Power [dB]')
        else:
            try:
                cax.xaxis.offsetText.set_visible(False)
                offset = cax.xaxis.get_major_formatter()
                # order_of_mag = str(offset.get_offset()).split('e')[-1]
                # This will fail if offset is empty
                order_of_mag = int(np.log10(float(offset.get_offset())))
                cax.set_xlabel(r'Power [cts / $10^{%d}$]' % order_of_mag)
            except ValueError:
                cax.set_xlabel('Power [counts]')


        # --------
        axTimeseries = plt.subplot(gs[4:16, -4:])
    #         axTimeseries = plt.axes(rect_timeseries)
        print('Plotting Timeseries')
        self.plot_time_series(tight=True, f_start=f_start, f_stop=f_stop, orientation='v')
        axTimeseries.yaxis.set_major_formatter(nullfmt)
    #        axTimeseries.xaxis.set_major_formatter(nullfmt)

        # --------
        # Could exclude since it takes much longer to run than the other plots.
        if kurtosis:
            axKurtosis = plt.subplot(gs[16:20, :10])
    #             axKurtosis = plt.axes(rect_kurtosis)
            print('Plotting Kurtosis')
            self.plot_kurtosis(f_start=f_start, f_stop=f_stop)


        # --------
        axHeader = plt.subplot(gs[4:16, :8])
    #         axHeader = plt.axes(rect_header)
        print('Plotting Header')
        # Generate nicer header
        telescopes = {0: 'Fake data',
                      1: 'Arecibo',
                      2: 'Ooty',
                      3: 'Nancay',
                      4: 'Parkes',
                      5: 'Jodrell',
                      6: 'GBT',
                      8: 'Effelsberg',
                      10: 'SRT',
                      64: 'MeerKAT',
                      65: 'KAT7'
                      }


    #     plot_header = "{:10} {:>20}\n".format("Telescope:", telescope)
        plot_header = ''
        for key in ['SOURCE_NAME',
                    'TELESCOPE_ID',
                    'SRC_RAJ', 
                    'SRC_DEJ',
                    'TSTART', 
                    'NCHANS',
                    'NBEAMS', 
                    'NIFS', 
                    'NBITS',
                    'FOFF',
                    'FCH1']:
            try:
                if key.lower() in self.header.keys():
                    value = self.header[key.lower()]
                else:
                    value = self.header[key.lower().encode()]

                if key == 'SRC_RAJ':
                    h, m, s = self.header[b'src_raj'].hms
                    value = '%02dh %02dm %06.3fs' % (h, m, s)
                elif key == 'SRC_DEJ':
                    d, m, s = self.header[b'src_dej'].dms
                    m = abs(m)
                    s = abs(s)
                    if d < 0:
                        str_d = '-%02d' % (-d)
                    else:
                        str_d = '+%02d' % d
                    value = '%sd %02dm %06.3fs' % (str_d, m, s)
                elif key == 'FOFF':
                    value = (value * 1e6 * u.Hz)
                    if np.abs(value) > 1e6 * u.Hz:
                        value = '%.06f MHz' % (value.to('MHz').value)
                    elif np.abs(value) > 1e3 * u.Hz:
                #             foff = "%.3f kHz" % (foff.to('kHz').value)
                        value = '%.03f kHz' % (value.to('kHz').value)
                    else:
                        value = '%d Hz'% (value.to('Hz').value)
                elif key == 'TSTART':
                    value = '%.02f' % value
                    key = 'TSTART (MJD)'
                elif key == 'FCH1':
                    value = "%6.6f MHz" % value
                elif key == 'TELESCOPE_ID':
                    value = telescopes[value]
                elif key == 'SOURCE_NAME':
                    try:
                        value = value.decode()
                    except AttributeError:
                        pass

                if key in ['SRC_DEJ', 'FOFF']:
                    plot_header += "{:<14} {:<18}\n".format('%s:' % key, str(value))
                else:
                    plot_header += "{:<15} {:<17}\n".format('%s:' % key, str(value))
            except KeyError:
                pass

        fig = plt.gcf()
        size = fig.get_size_inches()

        fontsize = int(size[0])
        plt.text(-0.1, 1, plot_header, fontname='monospace', fontsize=fontsize, ha='left', va='top', wrap=True)
        
        axHeader.axis('off')

        plt.tight_layout()

    def write_to_filterbank(self, filename_out):
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

    def calibrate_band_pass_N1(self):
        """ One way to calibrate the band pass is to take the median value
            for every frequency fine channel, and divide by it.
        """

        band_pass = np.median(self.data.squeeze(),axis=0)
        self.data = self.data/band_pass


def cmd_tool(args=None):
    """ Command line tool for plotting and viewing info on filterbank files """

    from argparse import ArgumentParser

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
