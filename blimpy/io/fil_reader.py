r'''Reader for Filterbank files (.fil)'''

import os

import numpy as np

from blimpy.io import sigproc
from blimpy.io.base_reader import Reader, logger, GIGA


class FilReader(Reader):
    """ This class handles .fil files.
    """

    def __init__(self, filename,f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True, max_load=None):
        """ Constructor.

        Args:
            filename (str): filename of blimpy file.
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            t_start (int): start time bin
            t_stop (int): stop time bin
            max_load (float): memory limit in gigabytes
        """
        super(FilReader, self).__init__()

        self.header_keywords_types = sigproc.header_keyword_types

        if filename and os.path.isfile(filename):
            self.filename = filename
            self.load_data = load_data
            self.header = self.read_header()
            self.file_size_bytes = os.path.getsize(self.filename)
            self.idx_data = sigproc.len_header(self.filename)
            self.n_channels_in_file  = self.header['nchans']
            self.n_beams_in_file = self.header['nifs'] #Placeholder for future development.
            self.n_pols_in_file = 1 #Placeholder for future development.
            self._n_bytes = int(self.header['nbits'] / 8)  #number of bytes per digit.
            self._d_type = self._setup_dtype()
            self._setup_n_ints_in_file()
            self.file_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header['foff'] < 0:
                self.f_end  = self.header['fch1']
                self.f_begin  = self.f_end + self.n_channels_in_file*self.header['foff']
            else:
                self.f_begin  = self.header['fch1']
                self.f_end  = self.f_begin + self.n_channels_in_file*self.header['foff']

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
            if max_load is not None and max_load > 0:
                self.max_data_array_size = max_load * GIGA

            if self.file_size_bytes > self.max_data_array_size:
                self.large_file = True
            else:
                self.large_file = False

            if self.load_data:
                if self.large_file:
                    if self.f_start or self.f_stop or self.t_start or self.t_stop:
                        if self.isheavy():
                            self.warn_memory("Selection", self._calc_selection_size())
                            self._init_empty_selection()
                        else:
                            self.read_data()
                    else:
                        self.warn_memory("File", self.file_size_bytes)
                        self._init_empty_selection()
                else:
                    self.read_data()
            else:
                logger.debug("Skipping loading data ...")
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
            self.warn_memory("Selection", self._calc_selection_size())
            self.data = np.array([0], dtype=self._d_type)
            return None

        #Convert input frequencies into what their corresponding channel number would be.
        self._setup_chans()
        #Update frequencies ranges from channel number.
        self._setup_freqs()

        n_chans = self.header['nchans']
        n_chans_selected = self.selection_shape[self.freq_axis]
        n_ifs   = self.header['nifs']

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
#                     if self.header['foff'] < 0:
#                         dd = dd[::-1]

                self.data[ii, jj] = dd

                f.seek(int(self._n_bytes  * (n_chans - self.chan_stop_idx)), 1)  # Seek to start of next block

        # Give the FD back to the O/S.
        f.close()

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

#         if self.header['foff'] < 0:
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
