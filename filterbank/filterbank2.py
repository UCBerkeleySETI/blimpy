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
import time
import struct
import numpy as np
from pprint import pprint

from astropy import units as u
from astropy.coordinates import Angle
import scipy.stats
from matplotlib.ticker import NullFormatter
import h5py

import file_wrapper as fw
from utils import db, lin, rebin, closest
import bitshuffle.h5

import pdb;# pdb.set_trace()

# Check if $DISPLAY is set (for handling plotting on remote machines with no X-forwarding)
if os.environ.has_key('DISPLAY'):
    import pylab as plt
else:
    import matplotlib
    matplotlib.use('Agg')
    import pylab as plt


#------
# Logging set up
import logging
logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    format = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    format = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'

logging.basicConfig(format=format,stream=stream,level = level_log)


###
# Config values
###

MAX_PLT_POINTS      = 65536                  # Max number of points in matplotlib plot
MAX_IMSHOW_POINTS   = (8192, 4096)           # Max number of points in imshow plot
MAX_HEADER_BLOCKS   = 100                    # Max size of header (in 512-byte blocks)
MAX_BLOB_MB         = 2   #256          # Max size of blob in MB

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

    keyword = str(keyword)

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
        elif keyword not in header_keyword_types.keys():
            pass
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

    def __init__(self, filename=None, f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True,header_dict=None, data_array=None):
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
            self.ext = filename.split(".")[-1].strip().lower()  #File extension
            self.container = fw.open_file(filename, f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop,load_data=load_data)
            self.header = self.container.header
            self.n_ints_in_file = self.container.n_ints_in_file
            self.__setup_time_axis()
            self.heavy =  self.container.heavy
            self.file_shape = self.container.file_shape
            self.file_size_bytes = self.container.file_size_bytes
            self.selection_shape = self.container.selection_shape
            self.n_channels_in_file = self.container.n_channels_in_file

            self.__load_data()

        elif header_dict is not None and data_array is not None:
            self.gen_from_header(header_dict, data_array)
        else:
            pass

    def __load_data(self):
        '''
        '''

        self.data = self.container.data
        self.freqs = self.container.freqs
        self.timestamps = self.container.timestamps

    def gen_from_header(self, header_dict, data_array, f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True):
        self.filename = ''
        self.header = header_dict
        self.data = data_array
        self.n_ints_in_file = 0

        self._setup_freqs()

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

    def __setup_time_axis(self,t_start=None, t_stop=None):
        '''  Setup time axis.
        '''

        # now check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if t_start:
            ii_start = t_start
        if t_stop:
            ii_stop = t_stop
        n_ints = ii_stop - ii_start

        ## Setup time axis
        t0 = self.header['tstart']
        t_delt = self.header['tsamp']
        self.timestamps = np.arange(0, n_ints) * t_delt / 24./60./60 + t0

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None):
        ''' Reads data selection if small enough.
        '''

        self.container.read_data(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)

        self.__load_data()

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
        if self.data is not None:
            print "%16s : %32s" % ("Data shape", self.file_shape)
        if self.freqs is not None:
            print "%16s : %32s" % ("Start freq (MHz)", self.freqs[0])
            print "%16s : %32s" % ("Stop freq (MHz)", self.freqs[-1])

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

    def calc_n_coarse_chan(self):
        ''' This makes an attempt to calculate the number of coarse channels in a given file.
            It assumes for now that a single coarse channel is 2.9296875 MHz
        '''

        # Could add a telescope based coarse channel bandwith, or other discriminative.
        # if telescope_id == 'GBT':
        # or actually as is currently
        # if self.header['telescope_id'] == 6:

        coarse_chan_bw = 2.9296875

        bandwith = abs(self.header['nchans']*self.header['foff'])
        n_coarse_chan = int(bandwith / coarse_chan_bw)

        return n_coarse_chan

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
        ''' Plot the time series.

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
        plt.text(0.05,.95,plot_header,ha='left', va='top', wrap=True)

        axHeader.set_axis_bgcolor('white')
        axHeader.xaxis.set_major_formatter(nullfmt)
        axHeader.yaxis.set_major_formatter(nullfmt)

    def write_to_fil(self, filename_out):
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
            It check the file size then decides how to write the file.

        Args:
            filename_out (str): Name of output file
        """

        if self.heavy:
            self.__write_to_hdf5_heavy(filename_out)
        else:
            self.__write_to_hdf5_light(filename_out)

    def __write_to_hdf5_heavy(self, filename_out, *args, **kwargs):
        """ Write data to HDF5 file.

        Args:
            filename_out (str): Name of output file
        """

        t0 = time.time()
        block_size = 0

        #Note that I make the intentional difference between a chunk and a blob here.
        chunk_dim = self.__get_chunk_dimentions()
        blob_dim = self.__get_blob_dimentions(chunk_dim)
        n_blobs = self.container.calc_n_blobs(blob_dim)

        with h5py.File(filename_out, 'w') as h5:

#            h5.attrs['CLASS'] = 'FILTERBANK'

            dset = h5.create_dataset('data',
                              shape=self.file_shape,
                              chunks=chunk_dim,
#                               compression=bitshuffle.h5.H5FILTER,
#                               compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
                              compression="lzf",
                              dtype=self.data.dtype)

            dset_mask = h5.create_dataset('mask',
                                     shape=self.file_shape,
                                     chunks=chunk_dim,
#                                      compression=bitshuffle.h5.H5FILTER,
#                                      compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
                              compression="lzf",
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

            if blob_dim[2] < self.n_channels_in_file:

                logger.info('Filling in with data over %i n_blobs...'% n_blobs)
                for ii in range(0, n_blobs):
                    logger.info('Reading %i of %i' % (ii + 1, n_blobs))

                    bob = self.container.read_blob(blob_dim,n_blob=ii)

                    # Reverse array if frequency axis is flipped
                    c_start = self.container.c_start() + ii*blob_dim[2]
                    t_start = self.container.t_start + (c_start/self.n_channels_in_file)*blob_dim[0]
                    t_stop = t_start + blob_dim[0]

                    if self.header['foff'] < 0:
                        c_start = self.n_channels_in_file - (c_start)%self.n_channels_in_file
                        c_stop = c_start - blob_dim[2]
                    else:
                        c_start = (c_start)%self.n_channels_in_file
                        c_stop = c_start + blob_dim[2]

                    logger.debug(t_start,t_stop,c_start,c_stop)

                    dset[t_start:t_stop,0,c_start:c_stop] = bob[:]

            else:

                logger.info('Filling in with data over %i n_blobs...'% n_blobs)
                for ii in range(0, n_blobs):
                    logger.info('Reading %i of %i' % (ii + 1, n_blobs))

                    bob = self.container.read_blob(blob_dim,n_blob=ii)
                    t_start = self.container.t_start + ii*blob_dim[0]
                    t_stop = min((ii+1)*blob_dim[0],self.n_ints_in_file)

                    dset[t_start:t_stop] = bob[:]

        t1 = time.time()
        logger.info('Conversion time: %2.2fsec' % (t1- t0))

    def __write_to_hdf5_light(self, filename_out, *args, **kwargs):
        """ Write data to HDF5 file in one go.

        Args:
            filename_out (str): Name of output file
        """


        with h5py.File(filename_out, 'w') as h5:

            dset = h5.create_dataset('data',
                              data=self.data,
                              compression='lzf')

            dset_mask = h5.create_dataset('mask',
                                     shape=self.file_shape,
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

    def __get_blob_dimentions(self,chunk_dim):
        ''' Sets the blob dimmentions, trying to read around 256 MiB at a time. This is assuming chunk is about 1 MiB.
        '''

        freq_axis_size = min(self.n_channels_in_file,chunk_dim[2]*MAX_BLOB_MB)
        time_axis_size = chunk_dim[0] * MAX_BLOB_MB * chunk_dim[2] / freq_axis_size

        blob_dim = (time_axis_size, 1, freq_axis_size)

        return blob_dim

    def __get_chunk_dimentions(self):
        ''' Sets the chunking dimmentions depending on the file type.
        '''

        if 'gpuspec.0000.' in self.filename:
            logger.info('Detecting high frequency resolution data.')
            chunk_dim = (1,1,1048576)
            return chunk_dim
        elif 'gpuspec.0001.' in self.filename:
            logger.info('Detecting high time resolution data.')
            chunk_dim = (512,1,2048)
            return chunk_dim
        elif 'gpuspec.0002.' in self.filename:
            logger.info('Detecting intermediate frequency and time resolution data.')
            chunk_dim = (10,1,65536)
#            chunk_dim = (1,1,65536/4)
            return chunk_dim
        else:
            logger.warning('File format not know. Will use autoblobing.')
            chunk_dim = True
            return chunk_dim

def cmd_tool(args=None):
    """ Command line tool for plotting and viewing info on filterbank files """

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Command line utility for reading and plotting filterbank files.")

    parser.add_argument('filename', type=str,
                        help='Name of file to read')
    parser.add_argument('-p', action='store',  default='a', dest='what_to_plot', type=str,
                        help='Show: "w" waterfall (freq vs. time) plot; "s" integrated spectrum plot, \
                             "a" for all available plots and information; and more.')
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
    parser.add_argument('-H', action='store_true', default=False, dest='to_hdf5',
                       help='Write file to hdf5 format.')
    parser.add_argument('-F', action='store_true', default=False, dest='to_fil',
                       help='Write file to .fil format.')
    parser.add_argument('-o', action='store', default=None, dest='filename_out', type=str,
                        help='Filename output (if not probided, the name will be the same but with apropiate extension).')

    parse_args = parser.parse_args()

    # Open filterbank data
    filename = parse_args.filename
    load_data = not parse_args.info_only
    info_only = parse_args.info_only
    filename_out = parse_args.filename_out

    # only load one integration if looking at spectrum
    wtp = parse_args.what_to_plot
    if not wtp or 's' in wtp:
        if parse_args.t_start == None:
            t_start = 0
        else:
            t_start = parse_args.t_start
        t_stop  = t_start + 1

        if parse_args.average:
            t_start = None
            t_stop  = None
    else:
        t_start = parse_args.t_start
        t_stop  = parse_args.t_stop

    fil = Filterbank(filename, f_start=parse_args.f_start, f_stop=parse_args.f_stop,t_start=parse_args.t_start, t_stop=parse_args.t_stop,load_data=load_data)
    fil.info()

    #Check the size of selection.

    if fil.heavy or parse_args.to_hdf5 or parse_args.to_fil:
        info_only = True

    # And if we want to plot data, then plot data.

    if not info_only:
        print ''

        # check start & stop frequencies make sense
        #try:
        #    if parse_args.f_start:
        #        print "Start freq: %2.2f" % parse_args.f_start
        #        assert parse_args.f_start >= fil.freqs[0] or np.isclose(parse_args.f_start, fil.freqs[0])
        #
        #    if parse_args.f_stop:
        #        print "Stop freq: %2.2f" % parse_args.f_stop
        #        assert parse_args.f_stop <= fil.freqs[-1] or np.isclose(parse_args.f_stop, fil.freqs[-1])
        #except AssertionError:
        #    print "Error: Start and stop frequencies must lie inside file's frequency range."
        #    print "i.e. between %2.2f-%2.2f MHz." % (fil.freqs[0], fil.freqs[-1])
        #    exit()

        if parse_args.blank_dc:
            print "Blanking DC bin"
            n_coarse_chan = fil.calc_n_coarse_chan()
            fil.blank_dc(n_coarse_chan)

        if "w" in parse_args.what_to_plot:
            plt.figure("waterfall", figsize=(8, 6))
            fil.plot_waterfall(f_start=parse_args.f_start, f_stop=parse_args.f_stop)
        elif "s" in parse_args.what_to_plot:
            plt.figure("Spectrum", figsize=(8, 6))
            fil.plot_spectrum(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all')
        elif "mm" in parse_args.what_to_plot:
            plt.figure("min max", figsize=(8, 6))
            fil.plot_spectrum_min_max(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all')
        elif "k" in parse_args.what_to_plot:
            plt.figure("kurtosis", figsize=(8, 6))
            fil.plot_kurtosis(f_start=parse_args.f_start, f_stop=parse_args.f_stop)
        elif "t" in parse_args.what_to_plot:
            plt.figure("Time Series", figsize=(8, 6))
            fil.plot_time_series(f_start=parse_args.f_start, f_stop=parse_args.f_stop)
        elif "a" in parse_args.what_to_plot:
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all')

        if parse_args.plt_filename != '':
            plt.savefig(parse_args.plt_filename)

        if not parse_args.save_only:
            if os.environ.has_key('DISPLAY'):
                plt.show()
            else:
                print "No $DISPLAY available."


    else:

        if parse_args.to_hdf5 and parse_args.to_fil:
            raise warning('Either provide to_hdf5 or to_fil, but not both.')

        if parse_args.to_hdf5:
            if not filename_out:
                filename_out = filename.replace('.fil','.h5')
            elif '.h5' not in filename_out:
                filename_out = filename_out.replace('.fil','')+'.h5'

            print 'Writing file : %s'%(filename_out)
            fil.write_to_hdf5(filename_out)
            print 'File written.'

        if parse_args.to_fil:
            if not filename_out:
                filename_out = filename.replace('.h5','.fil')
            elif '.fil' not in filename_out:
                filename_out = filename_out.replace('.h5','')+'.fil'

            print 'Writing file : %s'%(filename_out)
            fil.write_to_fil(filename_out)
            print 'File written.'


if __name__ == "__main__":
    cmd_tool()
