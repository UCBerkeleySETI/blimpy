import os

import h5py
import numpy as np
from astropy.coordinates import Angle

from blimpy.io.base_reader import Reader, logger, MAX_DATA_ARRAY_SIZE_UNIT


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
            self.h5 = h5py.File(self.filename, mode='r')
            self.read_header()
            self.file_size_bytes = os.path.getsize(self.filename)  # In bytes
            self.n_ints_in_file = self.h5["data"].shape[self.time_axis] #
            self.n_channels_in_file = self.h5["data"].shape[self.freq_axis] #
            self.n_beams_in_file = self.header['nifs'] #Placeholder for future development.
            self.n_pols_in_file = 1 #Placeholder for future development.
            self._n_bytes = int(self.header['nbits'] / 8)  #number of bytes per digit.
            self._d_type = self._setup_dtype()
            self.file_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header['foff'] < 0:
                self.f_end = self.header['fch1']
                self.f_begin = self.f_end + self.n_channels_in_file*self.header['foff']
            else:
                self.f_begin = self.header['fch1']
                self.f_end = self.f_begin + self.n_channels_in_file*self.header['foff']

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
            #if six.PY3:
            #    key = bytes(key, 'ascii')
            if isinstance(val, bytes):
                val = val.decode('ascii')
            if key == 'src_raj':
                self.header[key] = Angle(val, unit='hr')
            elif key == 'src_dej':
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

#         if self.header['foff'] < 0:
#             blob = blob[:,:,::-1]

        return blob