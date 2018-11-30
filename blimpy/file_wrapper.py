#!/usr/bin/env python
""" This modele handles file types.
"""

import os
import sys
import numpy as np
import h5py
import six

from astropy.coordinates import Angle

try:
    from . import sigproc
except:
    import sigproc

# import pdb;# pdb.set_trace()

import logging
logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    lformat = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    lformat = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'

logging.basicConfig(format=lformat, stream=stream, level=level_log)

# Max size of data array to load into memory (1GB in bytes)
MAX_DATA_ARRAY_SIZE_UNIT = 1024 * 1024 * 1024.0

class Reader(object):
    """ Basic reader object """

    def __init__(self):

        self.t_begin = 0
        self.t_end   = 0

    def _setup_selection_range(self, f_start=None, f_stop=None, t_start=None, t_stop=None, init=False):
        """Making sure the selection if time and frequency are within the file limits.

        Args:
            init (bool): If call during __init__
        """

        # This avoids resetting values
        if init is True:
            if t_start is None:
                t_start = self.t_begin
            if t_stop is None:
                t_stop = self.t_end
            if f_start is None:
                f_start = self.f_begin
            if f_stop is None:
                f_stop = self.f_end
        else:
            if f_start is None:
                f_start = self.f_start
            if f_stop is None:
                f_stop = self.f_stop
            if t_start is None:
                t_start = self.t_start
            if t_stop is None:
                t_stop = self.t_stop

        # By now, all values start/stop are populated.

        if t_stop >= 0 and t_start >= 0 and t_stop < t_start:
            t_stop, t_start = t_start,t_stop
            logger.warning('Given t_stop < t_start, assuming reversed values.')
        if f_stop and f_start and f_stop < f_start:
            f_stop, f_start = f_start,f_stop
            logger.warning('Given f_stop < f_start, assuming reversed values.')

        if t_start >= self.t_begin and t_start < self.t_end:
            self.t_start = int(t_start)
        else:
            if init is False or t_start != None:
                logger.warning('Setting t_start = %f, since t_start not given or not valid.'%self.t_begin)
            self.t_start = self.t_begin

        if t_stop <= self.t_end  and t_stop > self.t_begin:
            self.t_stop = int(t_stop)
        else:
            if init is False or t_stop:
                logger.warning('Setting t_stop = %f, since t_stop not given or not valid.'%self.t_end)
            self.t_stop = self.t_end

        if f_start >= self.f_begin and f_start < self.f_end:
            self.f_start = f_start
        else:
            if init is False or f_start:
                logger.warning('Setting f_start = %f, since f_start not given or not valid.'%self.f_begin)
            self.f_start = self.f_begin

        if f_stop <= self.f_end and f_stop > self.f_begin:
            self.f_stop = f_stop
        else:
            if init is False or f_stop:
                logger.warning('Setting f_stop = %f, since f_stop not given or not valid.'%self.f_end)
            self.f_stop = self.f_end

        # Now we have setup bounds, we can calculate shape of selection
        self.selection_shape = self._calc_selection_shape()

    def _init_empty_selection(self):
        """
        """

        self.data = np.array([0],dtype=self._d_type)

    def _setup_dtype(self):
        """Calculating dtype
        """

        #Set up the data type
        if self._n_bytes  == 4:
            return b'float32'
        elif self._n_bytes  == 2:
            return b'uint16'
        elif self._n_bytes  == 1:
            return b'uint8'
        else:
            logger.warning('Having trouble setting dtype, assuming float32.')
            return b'float32'

    def _calc_selection_size(self):
        """Calculate size of data of interest.
        """

        #Check to see how many integrations requested
        n_ints = self.t_stop - self.t_start
        #Check to see how many frequency channels requested
        n_chan = (self.f_stop - self.f_start) / abs(self.header[b'foff'])

        n_bytes  = self._n_bytes
        selection_size = int(n_ints*n_chan*n_bytes)

        return selection_size

    def _calc_selection_shape(self):
        """Calculate shape of data of interest.
        """

        #Check how many integrations requested
        n_ints = int(self.t_stop - self.t_start)
        #Check how many frequency channels requested
        n_chan = int(np.round((self.f_stop - self.f_start) / abs(self.header[b'foff'])))

        selection_shape = (n_ints,int(self.header[b'nifs']),n_chan)

        return selection_shape

    def _setup_chans(self):
        """Setup channel borders
        """

        if self.header[b'foff'] < 0:
            f0 = self.f_end
        else:
            f0 = self.f_begin

        i_start, i_stop = 0, self.n_channels_in_file
        if self.f_start:
            i_start = np.round((self.f_start - f0) / self.header[b'foff'])
        if self.f_stop:
            i_stop  = np.round((self.f_stop - f0)  / self.header[b'foff'])

        #calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx  = np.int(i_stop)

        if chan_stop_idx < chan_start_idx:
            chan_stop_idx, chan_start_idx = chan_start_idx,chan_stop_idx

        self.chan_start_idx =  chan_start_idx
        self.chan_stop_idx = chan_stop_idx

    def _setup_freqs(self):
        """Updating frequency borders from channel values
        """

        if self.header[b'foff'] > 0:
            self.f_start = self.f_begin + self.chan_start_idx*abs(self.header[b'foff'])
            self.f_stop = self.f_begin + self.chan_stop_idx*abs(self.header[b'foff'])
        else:
            self.f_start = self.f_end - self.chan_stop_idx*abs(self.header[b'foff'])
            self.f_stop = self.f_end - self.chan_start_idx*abs(self.header[b'foff'])

    def populate_timestamps(self,update_header=False):
        """  Populate time axis.
            IF update_header then only return tstart
        """

        #Check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if self.t_start:
            ii_start = self.t_start
        if self.t_stop:
            ii_stop = self.t_stop

        ## Setup time axis
        t0 = self.header[b'tstart']
        t_delt = self.header[b'tsamp']

        if update_header:
            timestamps = ii_start * t_delt / 24./60./60. + t0
        else:
            timestamps = np.arange(ii_start, ii_stop) * t_delt / 24./60./60. + t0

        return timestamps

    def populate_freqs(self):
        """
         Populate frequency axis
        """

        if self.header[b'foff'] < 0:
            f0 = self.f_end
        else:
            f0 = self.f_begin

        self._setup_chans()

        #create freq array
        i_vals = np.arange(self.chan_start_idx, self.chan_stop_idx)
        freqs = self.header[b'foff'] * i_vals + f0

        return freqs

    def calc_n_coarse_chan(self, chan_bw=None):
        """ This makes an attempt to calculate the number of coarse channels in a given file.

            Note:
                This is unlikely to work on non-Breakthrough Listen data, as a-priori knowledge of
                the digitizer system is required.

        """
        nchans = int(self.header[b'nchans'])
        # Do we have a file with enough channels that it has coarse channelization?
        if chan_bw is not None:
            bandwidth = abs(self.f_stop - self.f_start)
            n_coarse_chan = int(bandwidth / chan_bw)
            return n_coarse_chan
        elif nchans >= 2**20:
            # Does the common FFT length of 2^20 divide through without a remainder?
            # This should work for most GBT and all Parkes hires data
            if nchans % 2**20 == 0:
                n_coarse_chan = nchans // 2**20
                return n_coarse_chan
            # Early GBT data has non-2^N FFT length, check if it is GBT data
            elif self.header[b'telescope_id'] == 6:
                coarse_chan_bw = 2.9296875
                bandwidth = abs(self.f_stop - self.f_start)
                n_coarse_chan = int(bandwidth / coarse_chan_bw)
                return n_coarse_chan
            else:
                logger.warning("Couldn't figure out n_coarse_chan")
        elif self.header[b'telescope_id'] == 6 and nchans < 2**20:
            #For GBT non-hires data
            coarse_chan_bw = 2.9296875
            bandwidth = abs(self.f_stop - self.f_start)
            n_coarse_chan = int(bandwidth / coarse_chan_bw)
            return n_coarse_chan
        else:
            logger.warning("This function currently only works for hires BL Parkes or GBT data.")

    def calc_n_blobs(self, blob_dim):
        """ Given the blob dimensions, calculate how many fit in the data selection.
        """

        n_blobs = int(np.ceil(1.0 * np.prod(self.selection_shape) / np.prod(blob_dim)))

        return n_blobs

    def isheavy(self):
        """ Check if the current selection is too large.
        """

        selection_size_bytes = self._calc_selection_size()

        if selection_size_bytes > self.MAX_DATA_ARRAY_SIZE:
            return True
        else:
            return False


class H5Reader(Reader):
    """ This class handles .h5 files.
    """

    def __init__(self, filename, f_start=None, f_stop=None, t_start=None, t_stop=None, load_data=True, max_load=1.):
        """ Constructor.

        Args:
            filename (str): filename of blimpy file.
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            t_start (int): start time bin
            t_stop (int): stop time bin
        """
        super(H5Reader, self).__init__()

        if filename and os.path.isfile(filename) and h5py.is_hdf5(filename):

            #These values may be modified once code for multi_beam and multi_stokes observations are possible.
            self.freq_axis = 2
            self.time_axis = 0
            self.beam_axis = 1  # Place holder
            self.stokes_axis = 4  # Place holder

            self.filename = filename
            self.filestat = os.stat(filename)
            self.filesize = self.filestat.st_size/(1024.0**2)
            self.load_data = load_data
            self.h5 = h5py.File(self.filename)
            self.read_header()
            self.file_size_bytes = os.path.getsize(self.filename)  # In bytes
            self.n_ints_in_file = self.h5["data"].shape[self.time_axis] #
            self.n_channels_in_file = self.h5["data"].shape[self.freq_axis] #
            self.n_beams_in_file = self.header[b'nifs'] #Placeholder for future development.
            self.n_pols_in_file = 1 #Placeholder for future development.
            self._n_bytes = int(self.header[b'nbits'] / 8)  #number of bytes per digit.
            self._d_type = self._setup_dtype()
            self.file_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header[b'foff'] < 0:
                self.f_end = self.header[b'fch1']
                self.f_begin = self.f_end + self.n_channels_in_file*self.header[b'foff']
            else:
                self.f_begin = self.header[b'fch1']
                self.f_end = self.f_begin + self.n_channels_in_file*self.header[b'foff']

            self.t_begin = 0
            self.t_end = self.n_ints_in_file

            #Taking care all the frequencies are assigned correctly.
            self._setup_selection_range(f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop, init=True)
            #Convert input frequencies into what their corresponding channel number would be.
            self._setup_chans()
            #Update frequencies ranges from channel number.
            self._setup_freqs()

            #Applying data size limit to load.
            if max_load is not None:
                if max_load > 1.0:
                    logger.warning('Setting data limit > 1GB, please handle with care!')
                self.MAX_DATA_ARRAY_SIZE = max_load * MAX_DATA_ARRAY_SIZE_UNIT
            else:
                self.MAX_DATA_ARRAY_SIZE = MAX_DATA_ARRAY_SIZE_UNIT

            if self.file_size_bytes > self.MAX_DATA_ARRAY_SIZE:
                self.large_file = True
            else:
                self.large_file = False

            if self.load_data:
                if self.large_file:
                    #Only checking the selection, if the file is too large.
                    if self.f_start or self.f_stop or self.t_start or self.t_stop:
                        if self.isheavy():
                            logger.warning("Selection size of %.2f GB, exceeding our size limit %.2f GB. Instance created, header loaded, but data not loaded, please try another (t,v) selection." % (self._calc_selection_size() / (1024. ** 3), self.MAX_DATA_ARRAY_SIZE / (1024. ** 3)))
                            self._init_empty_selection()
                        else:
                            self.read_data()
                    else:
                        logger.warning("The file is of size %.2f GB, exceeding our size limit %.2f GB. Instance created, header loaded, but data not loaded. You could try another (t,v) selection."%(self.file_size_bytes/(1024.**3), self.MAX_DATA_ARRAY_SIZE/(1024.**3)))
                        self._init_empty_selection()
                else:
                    self.read_data()
            else:
                logger.info("Skipping loading data ...")
                self._init_empty_selection()
        else:
            raise IOError("Need a file to open, please give me one!")

    def read_header(self):
        """ Read header and return a Python dictionary of key:value pairs
        """

        self.header = {}

        for key, val in self.h5['data'].attrs.items():
            if six.PY3:
                key = bytes(key, 'ascii')
            if key == b'src_raj':
                self.header[key] = Angle(val, unit='hr')
            elif key == b'src_dej':
                self.header[key] = Angle(val, unit='deg')
            else:
                self.header[key] = val

        return self.header

    def _find_blob_start(self, blob_dim, n_blob):
        """Find first blob from selection.
        """

        #Convert input frequencies into what their corresponding channel number would be.
        self._setup_chans()

        #Check which is the blob time offset
        blob_time_start = self.t_start + blob_dim[self.time_axis]*n_blob

        #Check which is the blob frequency offset (in channels)
        blob_freq_start = self.chan_start_idx + (blob_dim[self.freq_axis]*n_blob)%self.selection_shape[self.freq_axis]

        blob_start = np.array([blob_time_start, 0, blob_freq_start])

        return blob_start

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None):
        """ Read data
        """

        self._setup_selection_range(f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop)

        #check if selection is small enough.
        if self.isheavy():
            logger.warning("Selection size of %.2f GB, exceeding our size limit %.2f GB. Instance created, header loaded, but data not loaded, please try another (t,v) selection." % (self._calc_selection_size() / (1024. ** 3), self.MAX_DATA_ARRAY_SIZE / (1024. ** 3)))
            self.data = np.array([0],dtype=self._d_type)
            return None

        #Convert input frequencies into what their corresponding channel number would be.
        self._setup_chans()
        #Update frequencies ranges from channel number.
        self._setup_freqs()

        self.data = self.h5["data"][self.t_start:self.t_stop,:,self.chan_start_idx:self.chan_stop_idx]

    def read_blob(self,blob_dim,n_blob=0):
        """Read blob from a selection.
        """

        n_blobs = self.calc_n_blobs(blob_dim)
        if n_blob > n_blobs or n_blob < 0:
            raise ValueError('Please provide correct n_blob value. Given %i, but max values is %i'%(n_blob,n_blobs))

        #This prevents issues when the last blob is smaller than the others in time
        if blob_dim[self.time_axis]*(n_blob+1) > self.selection_shape[self.time_axis]:
            updated_blob_dim = (self.selection_shape[self.time_axis] - blob_dim[self.time_axis]*n_blob, 1, blob_dim[self.freq_axis])
        else:
            updated_blob_dim = [int(i) for i in blob_dim]

        blob_start = self._find_blob_start(blob_dim, n_blob)
        blob_end = blob_start + np.array(updated_blob_dim)

        blob = self.h5["data"][int(blob_start[self.time_axis]):int(blob_end[self.time_axis]),
                               :,
                               int(blob_start[self.freq_axis]):int(blob_end[self.freq_axis])
                               ]

#         if self.header[b'foff'] < 0:
#             blob = blob[:,:,::-1]

        return blob

class FilReader(Reader):
    """ This class handles .fil files.
    """

    def __init__(self, filename,f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True, max_load=1.):
        """ Constructor.

        Args:
            filename (str): filename of blimpy file.
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            t_start (int): start time bin
            t_stop (int): stop time bin
        """
        super(FilReader, self).__init__()

        self.header_keywords_types = sigproc.header_keyword_types

        if filename and os.path.isfile(filename):
            self.filename = filename
            self.load_data = load_data
            self.header = self.read_header()
            self.file_size_bytes = os.path.getsize(self.filename)
            self.idx_data = sigproc.len_header(self.filename)
            self.n_channels_in_file  = self.header[b'nchans']
            self.n_beams_in_file = self.header[b'nifs'] #Placeholder for future development.
            self.n_pols_in_file = 1 #Placeholder for future development.
            self._n_bytes = int(self.header[b'nbits'] / 8)  #number of bytes per digit.
            self._d_type = self._setup_dtype()
            self._setup_n_ints_in_file()
            self.file_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header[b'foff'] < 0:
                self.f_end  = self.header[b'fch1']
                self.f_begin  = self.f_end + self.n_channels_in_file*self.header[b'foff']
            else:
                self.f_begin  = self.header[b'fch1']
                self.f_end  = self.f_begin + self.n_channels_in_file*self.header[b'foff']

            self.t_begin = 0
            self.t_end = self.n_ints_in_file

            #Taking care all the frequencies are assigned correctly.
            self._setup_selection_range(f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop, init=True)
            #Convert input frequencies into what their corresponding channel number would be.
            self._setup_chans()
            #Update frequencies ranges from channel number.
            self._setup_freqs()

            self.freq_axis = 2
            self.time_axis = 0
            self.beam_axis = 1  # Place holder

#EE ie.
#           spec = np.squeeze(fil_file.data)
            # set start of data, at real length of header  (future development.)
#            self.datastart=self.hdrraw.find('HEADER_END')+len('HEADER_END')+self.startsample*self.channels

            #Applying data size limit to load.
            if max_load is not None:
                if max_load > 1.0:
                    logger.warning('Setting data limit != 1GB, please handle with care!')
                self.MAX_DATA_ARRAY_SIZE = max_load * MAX_DATA_ARRAY_SIZE_UNIT
            else:
                self.MAX_DATA_ARRAY_SIZE = MAX_DATA_ARRAY_SIZE_UNIT

            if self.file_size_bytes > self.MAX_DATA_ARRAY_SIZE:
                self.large_file = True
            else:
                self.large_file = False

            if self.load_data:
                if self.large_file:
                    if self.f_start or self.f_stop or self.t_start or self.t_stop:
                        if self.isheavy():
                            logger.warning("Selection size of %.2f GB, exceeding our size limit %.2f GB. Instance created, header loaded, but data not loaded, please try another (t,v) selection." % (self._calc_selection_size() / (1024. ** 3), self.MAX_DATA_ARRAY_SIZE / (1024. ** 3)))
                            self._init_empty_selection()
                        else:
                            self.read_data()
                    else:
                        logger.warning("The file is of size %.2f GB, exceeding our size limit %.2f GB. Instance created, header loaded, but data not loaded. You could try another (t,v) selection."%(self.file_size_bytes/(1024.**3), self.MAX_DATA_ARRAY_SIZE/(1024.**3)))
                        self._init_empty_selection()
                else:
                    self.read_data()
            else:
                logger.info("Skipping loading data ...")
                self._init_empty_selection()
        else:
            raise IOError("Need a file to open, please give me one!")

    def _setup_n_ints_in_file(self):
        """ Calculate the number of integrations in the file. """
        self.n_ints_in_file = sigproc.calc_n_ints_in_file(self.filename)


    def read_header(self, return_idxs=False):
        """ Read blimpy header and return a Python dictionary of key:value pairs

        Args:
            filename (str): name of file to open

        Optional args:
            return_idxs (bool): Default False. If true, returns the file offset indexes
                                for values

        Returns:
            Python dict of key:value pairs, OR returns file offset indexes for values.

        """
        self.header = sigproc.read_header(self.filename, return_idxs=return_idxs)
        return self.header

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None):
        """ Read data.
        """

        self._setup_selection_range(f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop)

        #check if selection is small enough.
        if self.isheavy():
            logger.warning("Selection size of %.2f GB, exceeding our size limit %.2f GB. Instance created, "
                           "header loaded, but data not loaded, please try another (t,v) selection." % (self._calc_selection_size() / (1024. ** 3), self.MAX_DATA_ARRAY_SIZE / (1024. ** 3)))
            self.data = np.array([0],dtype=self._d_type)
            return None

        #Convert input frequencies into what their corresponding channel number would be.
        self._setup_chans()
        #Update frequencies ranges from channel number.
        self._setup_freqs()

        n_chans = self.header[b'nchans']
        n_chans_selected = self.selection_shape[self.freq_axis]
        n_ifs   = self.header[b'nifs']

        # Load binary data
        f = open(self.filename, 'rb')
        f.seek(int(self.idx_data))

        # now check to see how many integrations requested
        n_ints = self.t_stop - self.t_start

        # Seek to first integration
        f.seek(int(self.t_start * self._n_bytes  * n_ifs * n_chans), 1)

        #Loading  data
        self.data = np.zeros((n_ints, n_ifs, n_chans_selected), dtype=self._d_type)

        for ii in range(n_ints):
            for jj in range(n_ifs):
                f.seek(int(self._n_bytes  * self.chan_start_idx), 1) # 1 = from current location
                dd = np.fromfile(f, count=n_chans_selected, dtype=self._d_type)

                # Reverse array if frequency axis is flipped
#                     if self.header[b'foff'] < 0:
#                         dd = dd[::-1]

                self.data[ii, jj] = dd

                f.seek(int(self._n_bytes  * (n_chans - self.chan_stop_idx)), 1)  # Seek to start of next block

    def _find_blob_start(self):
        """Find first blob from selection.
        """

        # Convert input frequencies into what their corresponding channel number would be.
        self._setup_chans()

        # Check which is the blob time offset
        blob_time_start = self.t_start

        # Check which is the blob frequency offset (in channels)
        blob_freq_start = self.chan_start_idx

        blob_start = blob_time_start * self.n_channels_in_file + blob_freq_start

        return blob_start

    def read_blob(self,blob_dim,n_blob=0):
        """Read blob from a selection.
        """

        n_blobs = self.calc_n_blobs(blob_dim)
        if n_blob > n_blobs or n_blob < 0:
            raise ValueError('Please provide correct n_blob value. Given %i, but max values is %i'%(n_blob,n_blobs))

        # This prevents issues when the last blob is smaller than the others in time.
        if blob_dim[self.time_axis]*(n_blob+1) > self.selection_shape[self.time_axis]:
            updated_blob_dim = (int(self.selection_shape[self.time_axis] - blob_dim[self.time_axis]*n_blob), 1, int(blob_dim[self.freq_axis]))
        else:
            updated_blob_dim = [int(i) for i in blob_dim]

        blob_start = self._find_blob_start()
        blob = np.zeros(updated_blob_dim, dtype=self._d_type)

        # EE: For now; also assuming one polarization and one beam.

        # Assuming the blob will loop over the whole frequency range.
        if self.f_start == self.f_begin and self.f_stop == self.f_end:

            blob_flat_size = np.prod(blob_dim)
            updated_blob_flat_size = np.prod(updated_blob_dim)

            # Load binary data
            with open(self.filename, 'rb') as f:
                f.seek(int(self.idx_data + self._n_bytes  * (blob_start + n_blob*blob_flat_size)))
                dd = np.fromfile(f, count=updated_blob_flat_size, dtype=self._d_type)

            if dd.shape[0] == updated_blob_flat_size:
                blob = dd.reshape(updated_blob_dim)
            else:
                logger.info('DD shape != blob shape.')
                blob = dd.reshape((int(dd.shape[0]/blob_dim[self.freq_axis]),blob_dim[self.beam_axis],blob_dim[self.freq_axis]))
        else:

            for blobt in range(updated_blob_dim[self.time_axis]):

                #Load binary data
                with open(self.filename, 'rb') as f:
                    f.seek(int(self.idx_data + self._n_bytes * (blob_start + n_blob*blob_dim[self.time_axis]*self.n_channels_in_file + blobt*self.n_channels_in_file)))
                    dd = np.fromfile(f, count=blob_dim[self.freq_axis], dtype=self._d_type)

                blob[blobt] = dd

#         if self.header[b'foff'] < 0:
#             blob = blob[:,:,::-1]

        return blob

    def read_all(self,reverse=True):
        """ read all the data.
            If reverse=True the x axis is flipped.
        """
        raise NotImplementedError('To be implemented')

        # go to start of the data
        self.filfile.seek(int(self.datastart))
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
        self.filfile.seek(int(self.datastart+self.channels*rownumber*(int(self.nbits/8))))
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
        self.filfile.seek(int(self.datastart+self.channels*rownumber*(int(self.nbits/8))))
        # read data into 2-D numpy array
        data=np.fromfile(self.filfile,count=self.channels*n_rows,dtype=self.dtype).reshape(n_rows, self.channels)
        if reverse:
            data = data[:,::-1]
        return data


def open_file(filename, f_start=None, f_stop=None,t_start=None, t_stop=None,load_data=True,max_load=1.):
    """Open a HDF5 or filterbank file

    Returns instance of a Reader to read data from file.

    ================== ==================================================
    Filename extension File type
    ================== ==================================================
    h5, hdf5           HDF5 format
    fil                fil format
    *other*            Will raise NotImplementedError
    ================== ==================================================

    """
    if not os.path.isfile(filename):
        type(filename)
        print(filename)
        raise IOError("No such file or directory: " + filename)

    filename = os.path.expandvars(os.path.expanduser(filename))
    # Get file extension to determine type
    ext = filename.split(".")[-1].strip().lower()

    if six.PY3:
        ext = bytes(ext, 'ascii')

    if h5py.is_hdf5(filename):
        # Open HDF5 file
        return H5Reader(filename, f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop,
                        load_data=load_data, max_load=max_load)
    elif sigproc.is_filterbank(filename):
        # Open FIL file
        return FilReader(filename, f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop, load_data=load_data, max_load=max_load)
    else:
        raise NotImplementedError('Cannot open this type of file with Waterfall')
