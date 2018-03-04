#!/usr/bin/env python
""" This modele handles file types.
"""

import os
import sys
import struct
import numpy as np
import h5py

from astropy import units as u
from astropy.coordinates import Angle

import pdb;# pdb.set_trace()

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

MAX_DATA_ARRAY_SIZE = 1024 * 1024 * 1024.        # Max size of data array to load into memory (in bytes)


class  H5_reader(object):
    """ This class handles .h5 files.
    """

    def __init__(self, filename, f_start=None, f_stop=None, t_start=None, t_stop=None, load_data=True):
        """ Constructor.

        Args:
            filename (str): filename of blimpy file.
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            t_start (int): start time bin
            t_stop (int): stop time bin
        """

        if filename and os.path.isfile(filename) and h5py.is_hdf5(filename):
            self.filename = filename
            self.filestat = os.stat(filename)
            self.filesize = self.filestat.st_size/(1024.0**2)
            self.load_data = load_data
            self.h5 = h5py.File(self.filename)
            self.__read_header()
            self.file_size_bytes = os.path.getsize(self.filename)  # In bytes
            self.n_ints_in_file  = self.h5["data"].shape[0] #
            self.n_channels_in_file  = self.h5["data"].shape[2] #
            self.n_beams_in_file = self.header['nifs'] #Placeholder for future development.
            self.n_pols_in_file = 0 #Placeholder for future development.
            self.__n_bytes = self.header['nbits'] / 8  #number of bytes per digit.
            self.file_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header['foff'] < 0:
                self.f_end  = self.header['fch1']
                self.f_begin  = self.f_end + self.n_channels_in_file*self.header['foff']
            else:
                self.f_begin  = self.header['fch1']
                self.f_end  = self.f_begin + self.n_channels_in_file*self.header['foff']

            self.t_begin = 0
            self.t_end = self.n_ints_in_file

            self.__setup_selection_range(f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop,init=True)

            self.c_start = lambda: int(np.round((self.f_start - self.f_begin )/ abs(self.header['foff'])))
            self.c_stop = lambda: int(np.round((self.f_stop - self.f_begin )/ abs(self.header['foff'])))

            # These values will be modified once code for multi_beam and multi_stokes observations are possible.
            self.freq_axis = 2
            self.time_axis = 0
            self.beam_axis = 1  # Place holder
            self.stokes_axis = 4  # Place holder

            self.setup_time_axis()

            if self.file_size_bytes > MAX_DATA_ARRAY_SIZE:
                self.large_file = True
            else:
                self.large_file = False

            if self.load_data:
                if self.large_file:
                    #Only checking the selection, if the file is too large.
                    if self.f_start or self.f_stop or self.t_start or self.t_stop:
                        if self.isheavy():
                            logger.warning("Selection size of %.2f MB, exceeding our size limit %.2f MB. Instance created, header loaded, but data not loaded, please try another (t,v) selection."%(self.__calc_selection_size()/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                            self.data = np.array([0],dtype='float32')
                            self.freqs = np.array([0],dtype='float32')
                        else:
                            self.read_data()
                    else:
                        logger.warning("The file is of size %.2f MB, exceeding our size limit %.2f MB. Instance created, header loaded, but data not loaded. You could try another (t,v) selection."%(self.file_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                        self.data = np.array([0],dtype='float32')
                        self.freqs = np.array([0],dtype='float32')
                else:
                    self.read_data()
            else:
                print("Skipping loading data ...")
                self.data = np.array([0],dtype='float32')
                self.freqs = np.array([0],dtype='float32')

        else:
            raise IOError("Need a file to open, please give me one!")

    def __setup_selection_range(self, f_start=None, f_stop=None, t_start=None, t_stop=None,init=False):
        """Making sure the selection if time and frequency are within the file limits.

        Args:
            init (bool): If call during __init__
        """

        #This avoids reseting values
        if not init:
            if not f_start:
                f_start = self.f_start
            if not f_stop:
                f_stop = self.f_stop
            if not t_start:
                t_start = self.t_start
            if not t_stop:
                t_stop = self.t_stop

        if t_stop < t_start:
            t_stop, t_start = t_start,t_stop
            logger.warning('Given t_stop < t_start, assuming reversed values.')
        if f_stop < f_start:
            f_stop, f_start = f_start,f_stop
            logger.warning('Given f_stop < f_start, assuming reversed values.')

        if t_start and t_start >= self.t_begin and t_start < self.t_end:
            self.t_start = int(t_start)
        else:
            if not init:
                logger.warning('Setting t_start = t_begin.')
            self.t_start = self.t_begin

        if t_stop and t_stop <= self.t_end  and t_stop > self.t_begin:
            self.t_stop = int(t_stop)
        else:
            if not init:
                logger.warning('Setting t_stop = t_end.')
            self.t_stop = self.t_end

        if f_start and f_start >= self.f_begin and f_start < self.f_end:
            self.f_start = f_start
        else:
            if not init:
                logger.warning('Setting f_start = f_begin.')
            self.f_start = self.f_begin

        if f_stop and f_stop <= self.f_end and f_stop > self.f_begin:
            self.f_stop = f_stop
        else:
            if not init:
                logger.warning('Setting f_stop = f_end.')
            self.f_stop = self.f_end

        #calculate shape of selection
        self.selection_shape = self.__calc_selection_shape()

    def __read_header(self):
        """ Read header and return a Python dictionary of key:value pairs
        """

        self.header = {}

        for key, val in self.h5['data'].attrs.items():
            if key == 'src_raj':
                self.header[key] = Angle(val, unit='hr')
            elif key == 'src_dej':
                self.header[key] = Angle(val, unit='deg')
            else:
                self.header[key] = val

    def setup_time_axis(self,):
        """  Setup time axis.
        """

        #Check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if self.t_start:
            ii_start = self.t_start
        if self.t_stop:
            ii_stop = self.t_stop
        n_ints = ii_stop - ii_start

        ## Setup time axis
        t0 = self.header['tstart']
        t_delt = self.header['tsamp']
        self.timestamps = np.arange(0, n_ints) * t_delt / 24./60./60 + t0

    def __calc_selection_size(self):
        """Calculate size of data of interest.
        """

        #Check to see how many integrations requested
        n_ints = self.t_stop - self.t_start
        #Check to see how many frequency channels requested
        n_chan = (self.f_stop - self.f_start) / abs(self.header['foff'])

        n_bytes  = self.__n_bytes
        selection_size = int(n_ints*n_chan*n_bytes)

        return selection_size

    def isheavy(self):
        """ Check if the current selection is too large.
        """

        selection_size_bytes = self.__calc_selection_size()

        if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
            return True
        else:
            return False

    def __calc_selection_shape(self):
        """Calculate shape of data of interest.
        """

        #Check to see how many integrations requested
        n_ints = self.t_stop - self.t_start
        #Check to see how many frequency channels requested
        n_chan = int((self.f_stop - self.f_start) / abs(self.header['foff']))

        selection_shape = (n_ints,self.header['nifs'],n_chan)

        return selection_shape

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True):
        """ Read data
        """

        self.__setup_selection_range(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)

        #check if selection is small enough.
        if self.isheavy():
            logger.warning("Selection size of %.2f MB, exceeding our size limit %.2f MB. Instance created, header loaded, but data not loaded, please try another (t,v) selection."%(self.__calc_selection_size()/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            self.data = np.array([0],dtype='float32')
            self.freqs = np.array([0],dtype='float32')
            return None

        #convert input frequencies into what their corresponding index would be
        chan_start_idx, chan_stop_idx = self.__setup_freqs()

        self.data = self.h5["data"][self.t_start:self.t_stop,:,chan_start_idx:chan_stop_idx]

    def __setup_freqs(self):
        ## Setup frequency axis
        f0 = self.header['fch1']

        i_start, i_stop = 0, self.n_channels_in_file
        if self.f_start:
            i_start = (self.f_start - f0) / self.header['foff']
        if self.f_stop:
            i_stop  = (self.f_stop - f0)  / self.header['foff']

        #calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx  = np.int(i_stop)

        #create freq array
        if i_start < i_stop:
            i_vals = np.arange(chan_start_idx, chan_stop_idx)
        else:
            i_vals = np.arange(chan_stop_idx, chan_start_idx)

        self.freqs = self.header['foff'] * i_vals + f0

#         if self.header['foff'] < 0:
#             self.freqs = self.freqs[::-1]

        if chan_stop_idx < chan_start_idx:
            chan_stop_idx, chan_start_idx = chan_start_idx,chan_stop_idx

        return chan_start_idx, chan_stop_idx

    def calc_n_coarse_chan(self):
        """ This makes an attempt to calculate the number of coarse channels in a given file.
            It assumes for now that a single coarse channel is 2.9296875 MHz
        """

        # Could add a telescope based coarse channel bandwith, or other discriminative.
        # if telescope_id == 'GBT':
        # or actually as is currently
        # if self.header['telescope_id'] == 6:

        coarse_chan_bw = 2.9296875

        bandwith = abs(self.f_stop - self.f_start)
#        bandwith = abs(self.c_stop() - self.c_start())
        n_coarse_chan = int(bandwith / coarse_chan_bw)

        return max(n_coarse_chan, 1)

    def calc_n_blobs(self,blob_dim):
        """ Given the blob dimensions, calculate how many fit in the data selection.
        """

        n_blobs = int(np.ceil(self.__flat_array_dimmention(self.selection_shape)/float(self.__flat_array_dimmention(blob_dim))))

        return n_blobs

    def read_blob(self,blob_dim,n_blob=0):
        """Read blob from a selection.
        """

        n_blobs = self.calc_n_blobs(blob_dim)
        if n_blob > n_blobs or n_blob < 0:
            raise ValueError('Please provide correct n_blob value. Given %i, but max values is %i'%(n_blob,n_blobs))

        blob_start = self.__find_blob_start(blob_dim,n_blob)
        blob_end = blob_start + np.array(blob_dim)

#        blob = np.zeros(blob_dim,dtype='float32') #EE could remove.
        blob = self.h5["data"][blob_start[self.time_axis]:blob_end[self.time_axis],:,blob_start[self.freq_axis]:blob_end[self.freq_axis]]

#         if self.header['foff'] < 0:
#             blob = blob[:,:,::-1]

        return blob

    def __find_blob_start(self,blob_dim,n_blob):
        """Find first blob from selection.
        """

        #Check which is the blob time offset
        blob_time_start = self.t_start + blob_dim[self.time_axis]*n_blob

        #Check which is the blob frequency offset (in channels)
        blob_freq_start = self.c_start() + (blob_dim[self.freq_axis]*n_blob)%self.n_channels_in_file

        blob_start = np.array([blob_time_start,0,blob_freq_start])

        return blob_start

    def __flat_array_dimmention(self,array_dim):
        """Multiplies all the dimentions of an array.
        """

        array_flat_size = 1

        for a_dim in array_dim:
            array_flat_size*=a_dim

        return array_flat_size



class  FIL_reader(object):
    """ This class handles .fil files.
    """

    def __init__(self, filename,f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True):
        """ Constructor.

        Args:
            filename (str): filename of blimpy file.
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            t_start (int): start time bin
            t_stop (int): stop time bin
        """

        self.__set_header_keywords_types()

        if filename and os.path.isfile(filename):
            self.filename = filename
            self.load_data = load_data
            self.header = self.__read_header()
            self.file_size_bytes = os.path.getsize(self.filename)
            self.idx_data = self.__len_header()
            self.n_channels_in_file  = self.header['nchans']
            self.n_beams_in_file = self.header['nifs'] #Placeholder for future development.
            self.n_pols_in_file = 0 #Placeholder for future development.
            self.__n_bytes = self.header['nbits'] / 8  #number of bytes per digit.
            self.__get_n_ints_in_file()
            self.file_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header['foff'] < 0:
                self.f_end  = self.header['fch1']
                self.f_begin  = self.f_end + self.n_channels_in_file*self.header['foff']
            else:
                self.f_begin  = self.header['fch1']
                self.f_end  = self.f_begin + self.n_channels_in_file*self.header['foff']

            self.t_begin = 0
            self.t_end = self.n_ints_in_file

            self.__setup_selection_range(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop,init=True)

            self.c_start = lambda: int((self.f_start - self.f_begin )/ abs(self.header['foff']))
            self.c_stop = lambda: int((self.f_stop - self.f_begin )/ abs(self.header['foff']))

            self.setup_time_axis()

#EE ie.
#           spec = np.squeeze(fil_file.data)
            # set start of data, at real length of header  (future development.)
#            self.datastart=self.hdrraw.find('HEADER_END')+len('HEADER_END')+self.startsample*self.channels

            if self.file_size_bytes > MAX_DATA_ARRAY_SIZE:
                self.large_file = True
            else:
                self.large_file = False

            if self.load_data:
                if self.large_file:
                    if self.f_start or self.f_stop or self.t_start or self.t_stop:
                        if self.isheavy():
                            logger.warning("Selection size of %.2f MB, exceeding our size limit %.2f MB. Instance created, header loaded, but data not loaded, please try another (t,v) selection."%(self.__calc_selection_size()/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                            self.data = np.array([0],dtype='float32')
                            self.freqs = np.array([0],dtype='float32')
                        else:
                            self.read_data()
                    else:
                        logger.warning("The file is of size %.2f MB, exceeding our size limit %.2f MB. Instance created, header loaded, but data not loaded. You could try another (t,v) selection."%(self.file_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                        self.data = np.array([0],dtype='float32')
                        self.freqs = np.array([0],dtype='float32')
                else:
                    self.read_data()
            else:
                print("Skipping loading data ...")
                self.data = np.array([0],dtype='float32')
                self.freqs = np.array([0],dtype='float32')

        else:
            raise IOError("Need a file to open, please give me one!")

    def __setup_selection_range(self, f_start=None, f_stop=None, t_start=None, t_stop=None,init=False):
        """Making sure the selection if time and frequency are within the file limits.

            Args:
                init (bool): If call during __init__
        """

        #This avoids reseting values
        if not init:
            if not f_start:
                f_start = self.f_start
            if not f_stop:
                f_stop = self.f_stop
            if not t_start:
                t_start = self.t_start
            if not t_stop:
                t_stop = self.t_stop

        if t_stop < t_start:
            t_stop, t_start = t_start,t_stop
            logger.warning('Given t_stop < t_start, assuming reversed values.')
        if f_stop < f_start:
            f_stop, f_start = f_start,f_stop
            logger.warning('Given f_stop < f_start, assuming reversed values.')

        if t_start and t_start >= self.t_begin and t_start < self.t_end:
            self.t_start = int(t_start)
        else:
            if not init:
                logger.warning('Setting t_start = t_begin.')
            self.t_start = self.t_begin

        if t_stop and t_stop <= self.t_end  and t_stop > self.t_begin:
            self.t_stop = int(t_stop)
        else:
            if not init:
                logger.warning('Setting t_stop = t_end.')
            self.t_stop = self.t_end

        if f_start and f_start >= self.f_begin and f_start < self.f_end:
            self.f_start = f_start
        else:
            if not init:
                logger.warning('Setting f_start = f_begin.')
            self.f_start = self.f_begin

        if f_stop and f_stop <= self.f_end and f_stop > self.f_begin:
            self.f_stop = f_stop
        else:
            if not init:
                logger.warning('Setting f_stop = f_end.')
            self.f_stop = self.f_end

        #calculate shape of selection
        self.selection_shape = self.__calc_selection_shape()

    def __calc_selection_size(self):
        """Calculate size of data of interest.
        """

        #Check to see how many integrations requested
        n_ints = self.t_stop - self.t_start
        #Check to see how many frequency channels requested
        n_chan = (self.f_stop - self.f_start) / abs(self.header['foff'])

        n_bytes  = self.__n_bytes
        selection_size = int(n_ints*n_chan*n_bytes)

        return selection_size

    def isheavy(self):
        """ Check if the current selection is too large.
        """

        selection_size_bytes = self.__calc_selection_size()

        if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
            return True
        else:
            return False

    def __calc_selection_shape(self):
        """Calculate shape of data of interest.
        """

        #Check how many integrations were requested
        n_ints = self.t_stop - self.t_start
        #Check how many frequency channels were requested
        n_chan = int((self.f_stop - self.f_start) / abs(self.header['foff']))

        selection_shape = (n_ints,self.header['nifs'],n_chan)

        return selection_shape

    def __set_header_keywords_types(self):
        self.__header_keyword_types = {
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

    def __get_n_ints_in_file(self):

        n_bytes  = self.__n_bytes
        n_chans = self.n_channels_in_file
        n_ifs   = self.n_beams_in_file

        n_bytes_data = self.file_size_bytes - self.idx_data
        self.n_ints_in_file = n_bytes_data / (n_bytes * n_chans * n_ifs)

    def __len_header(self):
        """ Return the length of the blimpy header, in bytes

        Args:
            filename (str): name of file to open

        Returns:
            idx_end (int): length of header, in bytes
        """
        with  open(self.filename, 'rb') as f:
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

    def __read_next_header_keyword(self,fh):
        """

        Args:
            fh (file): file handler

        Returns:
        """
        n_bytes = np.fromstring(fh.read(4), dtype='uint32')[0]

        if n_bytes > 255:
            n_bytes = 16

        keyword = fh.read(n_bytes)

        #print keyword

        if keyword == 'HEADER_START' or keyword == 'HEADER_END':
            return keyword, 0, fh.tell()
        else:
            dtype = self.__header_keyword_types[keyword]
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
                val = self.__fil_double_to_angle(val)
                if keyword == 'src_raj':
                    val = Angle(val, unit=u.hour)
                else:
                    val = Angle(val, unit=u.deg)
            return keyword, val, idx

    def __read_header(self, return_idxs=False):
        """ Read blimpy header and return a Python dictionary of key:value pairs

        Args:
            filename (str): name of file to open

        Optional args:
            return_idxs (bool): Default False. If true, returns the file offset indexes
                                for values

        returns

        """
        with open(self.filename, 'rb') as fh:
            header_dict = {}
            header_idxs = {}

            # Check this is a blimpy file
            keyword, value, idx = self.__read_next_header_keyword(fh)

            try:
                assert keyword == 'HEADER_START'
            except AssertionError:
                raise RuntimeError("Not a valid blimpy file.")

            while True:
                keyword, value, idx = self.__read_next_header_keyword(fh)
                if keyword == 'HEADER_END':
                    break
                else:
                    header_dict[keyword] = value
                    header_idxs[keyword] = idx

        if return_idxs:
            return header_idxs
        else:
            return header_dict

    def __fil_double_to_angle(self,angle):
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

    def __setup_freqs(self,):
        ## Setup frequency axis
        f0 = self.header['fch1']

        i_start, i_stop = 0, self.n_channels_in_file
        if self.f_start:
            i_start = (self.f_start - f0) / self.header['foff']
        if self.f_stop:
            i_stop  = (self.f_stop - f0)  / self.header['foff']

        #calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx  = np.int(i_stop)

        #create freq array
        if i_start < i_stop:
            i_vals = np.arange(chan_start_idx, chan_stop_idx)
        else:
            i_vals = np.arange(chan_stop_idx, chan_start_idx)

        self.freqs = self.header['foff'] * i_vals + f0

#         if self.header['foff'] < 0:
#             self.freqs = self.freqs[::-1]

        return chan_start_idx, chan_stop_idx

    def setup_time_axis(self):
        """  Setup time axis.
        """

        #Check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if self.t_start:
            ii_start = self.t_start
        if self.t_stop:
            ii_stop = self.t_stop
        n_ints = ii_stop - ii_start

        ## Setup time axis
        t0 = self.header['tstart']
        t_delt = self.header['tsamp']
        self.timestamps = np.arange(0, n_ints) * t_delt / 24./60./60 + t0

    def info(self):
        """ Print header information and other derived information. """

        if not self.freqs:
            chan_start_idx, chan_stop_idx = self.__setup_freqs()

        for key, val in self.header.items():
            if key == 'src_raj':
                val = val.to_string(unit=u.hour, sep=':')
            if key == 'src_dej':
                val = val.to_string(unit=u.deg, sep=':')
            print("%16s : %32s" % (key, val))

        print("\n%16s : %32s" % ("Num ints in file", self.n_ints_in_file))
#        print "%16s : %32s" % ("Data shape", self.data.shape)
        print("%16s : %32s" % ("Start freq (MHz)", self.freqs[0]))
        print("%16s : %32s" % ("Stop freq (MHz)", self.freqs[-1]))

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True):
        """ Read data.
        """

        self.__setup_selection_range(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)

        #check if selection is small enough.
        if self.isheavy():
            logger.warning("Selection size of %.2f MB, exceeding our size limit %.2f MB. Instance created, header loaded, but data not loaded, please try another (t,v) selection."%(self.__calc_selection_size()/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            self.data = np.array([0],dtype='float32')
            self.freqs = np.array([0],dtype='float32')
            return None

        #convert input frequencies into what their corresponding index would be
        chan_start_idx, chan_stop_idx = self.__setup_freqs()

        n_chans = self.header['nchans']
        n_chans_selected = self.freqs.shape[0]
        n_ifs   = self.header['nifs']

        # Load binary data
        f = open(self.filename, 'rb')
        f.seek(self.idx_data)

        # now check to see how many integrations requested
        n_ints = self.t_stop - self.t_start

        # Seek to first integration
        f.seek(self.t_start * self.__n_bytes  * n_ifs * n_chans, 1)

        # Set up indexes used in file read (taken out of loop for speed)
        i0 = np.min((chan_start_idx, chan_stop_idx))
        i1 = np.max((chan_start_idx, chan_stop_idx))

        #Set up the data type (taken out of loop for speed)
        if self.__n_bytes  == 4:
            dd_type = 'float32'
        elif self.__n_bytes  == 2:
            dd_type = 'int16'
        elif self.__n_bytes  == 1:
            dd_type = 'int8'

        #EE Could add reading all in one go if file is small...then reshape. Unless actually reading a subsection.
        if load_data:
            self.data = np.zeros((n_ints, n_ifs, n_chans_selected), dtype='float32')

            for ii in range(n_ints):
                for jj in range(n_ifs):
                    f.seek(self.__n_bytes  * i0, 1) # 1 = from current location
                    dd = np.fromfile(f, count=n_chans_selected, dtype=dd_type)

                    # Reverse array if frequency axis is flipped
#                     if self.header['foff'] < 0:
#                         dd = dd[::-1]

                    self.data[ii, jj] = dd

                    f.seek(self.__n_bytes  * (n_chans - i1), 1)  # Seek to start of next block
        else:
            print("Skipping data load...")
            self.data = np.array([0])

        ## Setup time axis
        t0 = self.header['tstart']
        t_delt = self.header['tsamp']
        self.timestamps = np.arange(0, n_ints) * t_delt / 24./60./60 + t0

#    def __read_blob(self,blob_dim,n_blob=0):
    def read_blob(self,blob_dim,n_blob=0):
        """Read blob from a selection.
        """

        n_blobs = self.calc_n_blobs(blob_dim)
        if n_blob > n_blobs or n_blob < 0:
            raise ValueError('Please provide correct n_blob value. Given %i, but max values is %i'%(n_blob,n_blobs))
        blob_start = self.__find_blob_start(blob_dim)
        blob = np.zeros(blob_dim,dtype='float32')

        # Assuming the blob will be either one dimensional or loops over the whole frequency range.
        #EE: For now; also assuming one polarization and one beam.
        blob_flat_size = self.__flat_array_dimmention(blob_dim)

        # Load binary data
        f = open(self.filename, 'rb')
        f.seek(self.idx_data + self.__n_bytes  * (blob_start + n_blob*blob_flat_size))

        #Set up the data type (taken out of loop for speed)
        if self.__n_bytes  == 4:
            dd_type = 'float32'
        elif self.__n_bytes  == 2:
            dd_type = 'int16'
        elif self.__n_bytes  == 1:
            dd_type = 'int8'

        dd = np.fromfile(f, count=blob_flat_size, dtype=dd_type)

        if dd.shape[0] == blob_flat_size:
            blob = dd.reshape(blob_dim)
        else:
            logger.debug('DD shape != blob shape.')
            blob = dd.reshape((dd.shape[0]/blob_dim[2],blob_dim[1],blob_dim[2]))

#         if self.header['foff'] < 0:
#             blob = blob[:,:,::-1]

        return blob

    def __find_blob_start(self,blob_dim):
        """Find first blob from selection.
        """

        #Check which is the blob time offset
        blob_time_start = self.t_start

        #Check which is the blob frequency offset (in channels)
        blob_freq_start = self.c_start()

        blob_start = blob_time_start*self.n_channels_in_file + blob_freq_start

        return blob_start

#    def __calc_n_blobs(self,blob_dim):
    def calc_n_blobs(self,blob_dim):
        """ Given the blob dimensions, calculate how many fit in the data selection.
        """

        n_blobs = int(np.ceil(self.__flat_array_dimmention(self.selection_shape)/float(self.__flat_array_dimmention(blob_dim))))

        return n_blobs

    def __flat_array_dimmention(self,array_dim):
        """Multiplies all the dimentions of an array.
        """

        array_flat_size = 1

        for a_dim in array_dim:
            array_flat_size*=a_dim

        return array_flat_size

    def calc_n_coarse_chan(self):
        """ This makes an attempt to calculate the number of coarse channels in a given file.
            It assumes for now that a single coarse channel is 2.9296875 MHz
        """

        # Could add a telescope based coarse channel bandwith, or other discriminative.
        # if telescope_id == 'GBT':
        # or actually as is currently
        # if self.header['telescope_id'] == 6:

        coarse_chan_bw = 2.9296875

        bandwith = abs(self.f_stop - self.f_start)
        n_coarse_chan = int(bandwith / coarse_chan_bw)

        return n_coarse_chan

    def read_all(self,reverse=True):
        """ read all the data.
            If reverse=True the x axis is flipped.
        """
        raise NotImplementedError('To be implemented')

        # go to start of the data
        self.filfile.seek(self.datastart)
        # read data into 2-D numpy array
#        data=np.fromfile(self.filfile,dtype=self.dtype).reshape(self.channels,self.blocksize,order='F')
        data=np.fromfile(self.filfile,dtype=self.dtype).reshape(self.blocksize, self.channels)
        if reverse:
            data = data[:,::-1]
        return data

    def read_row(self,rownumber,reverse=True):
        """ Read a block of data. The number of samples per row is set in self.channels
            If reverse=True the x axis is flipped.
        """
        raise NotImplementedError('To be implemented')

        # go to start of the row
        self.filfile.seek(self.datastart+self.channels*rownumber*(self.nbits/8))
        # read data into 2-D numpy array
        data=np.fromfile(self.filfile,count=self.channels,dtype=self.dtype).reshape(1, self.channels)
        if reverse:
            data = data[:,::-1]
        return data

    def read_rows(self,rownumber,n_rows,reverse=True):
        """ Read a block of data. The number of samples per row is set in self.channels
            If reverse=True the x axis is flipped.
        """
        raise NotImplementedError('To be implemented')

        # go to start of the row
        self.filfile.seek(self.datastart+self.channels*rownumber*(self.nbits/8))
        # read data into 2-D numpy array
        data=np.fromfile(self.filfile,count=self.channels*n_rows,dtype=self.dtype).reshape(n_rows, self.channels)
        if reverse:
            data = data[:,::-1]
        return data


def open_file(filename, f_start=None, f_stop=None,t_start=None, t_stop=None,load_data=True):
    """Open a supported file type or fall back to Python built in open function.

    ================== ==================================================
    Filename extension File type
    ================== ==================================================
    h5                 HDF5 format
    fil                fil format
    *other*            Open with regular python :func:`open` function.
    ================== ==================================================

    """
    if not os.path.isfile(filename):
        type(filename)
        print(filename)
        raise IOError("No such directory: " + filename)

    filename = os.path.expandvars(os.path.expanduser(filename))
    # Get file extension to determine type
    ext = filename.split(".")[-1].strip().lower()

    if ext == 'h5':
        # Open HDF5 file
        return H5_reader(filename,f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop,load_data=load_data)
    elif ext == 'fil':
        # Open FIL file
        return FIL_reader(filename,f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop,load_data=load_data)
    else:
        # Fall back to regular Python `open` function
        return open(filename, *args, **kwargs)

