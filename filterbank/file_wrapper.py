#!/usr/bin/env python
''' This modele handles file types.
'''

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

#MAX_PLT_POINTS      = 65536                     # Max number of points in matplotlib plot
#MAX_IMSHOW_POINTS   = (8192, 4096)              # Max number of points in imshow plot
MAX_DATA_ARRAY_SIZE = 1024 * 1024 * 1024.        # Max size of data array to load into memory (in bytes)
#MAX_HEADER_BLOCKS   = 100                       # Max size of header (in 512-byte blocks)


class  H5_reader(object):
    ''' This class handles .h5 files.
    '''

    #EE check freq axis.


    def __init__(self, filename,f_start=None, f_stop=None,t_start=None, t_stop=None):
        """ Constructor.

        Args:
            filename (str): filename of filterbank file.
        """

        if filename and os.path.isfile(filename) and h5py.is_hdf5(filename):
            self.filename = filename
            self.data = None
            self.freqs = None
            self.h5 = h5py.File(self.filename)
            self.__read_header()
            self.file_size_bytes = os.path.getsize(self.filename)  # In bytes
            self.__setup_time_axis()
            self.n_ints_in_file  = self.h5["data"].shape[0] #
            self.n_channels_in_file  = self.h5["data"].shape[2] #
            self.n_beams_in_file = self.header['nifs'] #Placeholder for future development.
            self.n_pols_in_file = 0 #Placeholder for future development.
            self.__n_bytes = self.header['nbits'] / 8.  #number of bytes per digit.
            self.data_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header['foff'] < 0 :
                self.f_end  = self.header['foff']
                self.f_begin  = self.f_end + self.n_channels_in_file*self.header['foff']
            else:
                self.f_beging  = self.header['foff']
                self.f_end  = self.f_beging + self.n_channels_in_file*self.header['foff']

            self.t_start = t_start
            self.t_stop = t_stop
            self.f_start = f_start
            self.f_stop = f_stop
            self.c_start = (self.f_start - self.f_begin )/ self.header['foff']
            self.c_stop = (self.f_stop - self.f_begin )/ self.header['foff']

            if self.file_size_bytes > MAX_DATA_ARRAY_SIZE:
                self.heavy = True
                if self.f_start or self.f_stop or self.t_start or self.t_stop:
                    selection_size_bytes = self.__calc_selection_size()
                    if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
                        logger.warning("Selection size of %f MB, exceeding our size limit %f MB. Data not loaded, please try another (t,v) selection."%(selection_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                    else:
                        self.read_data()
                else:
                    logger.warning("The file is of size %f MB, exceeding our size limit %f MB. Data not loaded."%(self.file_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            else:
                self.heavy = False
                self.read_data()

        else:
            raise IOError("Need a file to open, please give me one!")

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

    def __setup_time_axis(self,):
        '''  Setup time axis.
        '''

        # now check to see how many integrations requested
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
        '''Calculate size of data of interest.
        '''

        #Check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if self.t_start:
            ii_start = self.t_start
        if self.t_stop:
            ii_stop = self.t_stop
        n_ints = ii_stop - ii_start

        #Check to see how many frequency channels requested
        jj_start, jj_stop = self.header['fch1']+self.header['foff']*self.header['nchans'], self.header['fch1']
        if self.f_start:
            jj_start = self.f_start
        if self.f_stop:
            jj_stop = self.f_stop
        n_chan = (jj_stop - jj_start) / abs(self.header['foff'])


        n_bytes  = self.__n_bytes

        selection_size = n_ints*n_chan*n_bytes

        return selection_size

    def read_data(self, load_data=True):
        ''' Read data
        '''


        #check if selection is small enough.
        selection_size_bytes = self.__calc_selection_size()
        if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
            logger.warning("Selection size of %f MB, exceeding our size limit %f MB. Data not loaded, please try another (t,v) selection."%(selection_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            return None


        self.data = self.h5["data"][self.t_start:self.t_stop,self.c_start:self.c_stop]
        self.__setup_freqs()

    def __setup_freqs(self):
        ## Setup frequency axis
        f0 = self.header['fch1']
        f_delt = self.header['foff']

        i_start, i_stop = 0, self.header['nchans']
        if self.f_start:
            i_start = (self.f_start - f0) / f_delt
        if self.f_stop:
            i_stop  = (self.f_stop - f0)  / f_delt

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


class  FIL_reader(object):
    ''' This class handles .fil files.
    '''

    def __init__(self, filename,f_start=None, f_stop=None,t_start=None, t_stop=None):
        """ Constructor.

        Args:
            filename (str): filename of filterbank file.
        """

        self.__set_header_keywords_types()

        if filename and os.path.isfile(filename):
            self.filename = filename
            self.data = None
            self.freqs = None
            self.header = self.__read_header()
            self.file_size_bytes = os.path.getsize(self.filename)
            self.idx_data = self.__len_header()
            self.__setup_time_axis()
            self.n_channels_in_file  = self.header['nchans']
            self.n_beams_in_file = self.header['nifs'] #Placeholder for future development.
            self.n_pols_in_file = 0 #Placeholder for future development.
            self.__n_bytes = self.header['nbits'] / 8.  #number of bytes per digit.
            self.__get_n_ints_in_file()
            self.data_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            if self.header['foff'] < 0 :
                self.f_end  = self.header['foff']
                self.f_begin  = self.f_end + self.n_channels_in_file*self.header['foff']
            else:
                self.f_beging  = self.header['foff']
                self.f_end  = self.f_beging + self.n_channels_in_file*self.header['foff']

            self.t_start = t_start
            self.t_stop = t_stop
            self.f_start = f_start
            self.f_stop = f_stop
            self.c_start = (self.f_start - self.f_begin )/ self.header['foff']
            self.c_stop = (self.f_stop - self.f_begin )/ self.header['foff']

#EE ie.
#           spec = np.squeeze(fil_file.data)
            # set start of data, at real length of header  (future development.)
#            self.datastart=self.hdrraw.find('HEADER_END')+len('HEADER_END')+self.startsample*self.channels

            if self.file_size_bytes > MAX_DATA_ARRAY_SIZE:
                self.heavy = True
                if self.f_start or self.f_stop or self.t_start or self.t_stop:
                    selection_size_bytes = self.__calc_selection_size()
                    if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
                        logger.warning("Selection size of %f MB, exceeding our size limit %f MB. Data not loaded, please try another (t,v) selection."%(selection_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                    else:
                        self.read_data()
                else:
                    logger.warning("The file is of size %f MB, exceeding our size limit %f MB. Data not loaded."%(self.file_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            else:
                self.heavy = False
                self.read_data()

        else:
            raise IOError("Need a file to open, please give me one!")

    def __calc_selection_size(self):
        '''Calculate size of data of interest.
        '''

        #Check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if self.t_start:
            ii_start = self.t_start
        if self.t_stop:
            ii_stop = self.t_stop
        n_ints = ii_stop - ii_start

        #Check to see how many frequency channels requested
        jj_start, jj_stop = self.header['fch1']+self.header['foff']*self.header['nchans'], self.header['fch1']
        if self.f_start:
            jj_start = self.f_start
        if self.f_stop:
            jj_stop = self.f_stop
        n_chan = (jj_stop - jj_start) / abs(self.header['foff'])


        n_bytes  = self.__n_bytes

        selection_size = n_ints*n_chan*n_bytes

        return selection_size

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
        """ Return the length of the filterbank header, in bytes

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
        """ Read filterbank header and return a Python dictionary of key:value pairs

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

            # Check this is a filterbank file
            keyword, value, idx = self.__read_next_header_keyword(fh)

            try:
                assert keyword == 'HEADER_START'
            except AssertionError:
                raise RuntimeError("Not a valid filterbank file.")

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
        f_delt = self.header['foff']

        i_start, i_stop = 0, self.header['nchans']
        if self.f_start:
            i_start = (self.f_start - f0) / f_delt
        if self.f_stop:
            i_stop  = (self.f_stop - f0)  / f_delt

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

    def __setup_time_axis(self):
        '''  Setup time axis.
        '''

        # now check to see how many integrations requested
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
            i_start, i_stop, chan_start_idx, chan_stop_idx = self.__setup_freqs()

        for key, val in self.header.items():
            if key == 'src_raj':
                val = val.to_string(unit=u.hour, sep=':')
            if key == 'src_dej':
                val = val.to_string(unit=u.deg, sep=':')
            print "%16s : %32s" % (key, val)

        print "\n%16s : %32s" % ("Num ints in file", self.n_ints_in_file)
#        print "%16s : %32s" % ("Data shape", self.data.shape)
        print "%16s : %32s" % ("Start freq (MHz)", self.freqs[0])
        print "%16s : %32s" % ("Stop freq (MHz)", self.freqs[-1])

    def read_data(self, load_data=True):
        ''' Read data.
        '''

        #check if selection is small enough.
        selection_size_bytes = self.__calc_selection_size()
        if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
            logger.warning("Selection size of %f MB, exceeding our size limit %f MB. Data not loaded, please try another (t,v) selection."%(selection_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            return None
#            load_data = False


        ## Setup frequency axis
        f0 = self.header['fch1']
        f_delt = self.header['foff']

        # keep this seperate!
        # file_freq_mapping =  np.arange(0, self.header['nchans'], 1, dtype='float64') * f_delt + f0

        #convert input frequencies into what their corresponding index would be

        i_start, i_stop, chan_start_idx, chan_stop_idx = self.__setup_freqs()

        n_bytes  = self.__n_bytes
        n_chans = self.header['nchans']
        n_chans_selected = self.freqs.shape[0]
        n_ifs   = self.header['nifs']

        # Load binary data
        f = open(self.filename, 'rb')
        f.seek(self.idx_data)

        # now check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if self.t_start:
            ii_start = self.t_start
        if self.t_stop:
            ii_stop = self.t_stop
        n_ints = ii_stop - ii_start

        # Seek to first integration
        f.seek(ii_start * n_bytes * n_ifs * n_chans, 1)

        # Set up indexes used in file read (taken out of loop for speed)
        i0 = np.min((chan_start_idx, chan_stop_idx))
        i1 = np.max((chan_start_idx, chan_stop_idx))

        #Set up the data type (taken out of loop for speed)
        if n_bytes == 4:
            dd_type = 'float32'
        elif n_bytes == 2:
            dd_type = 'int16'
        elif n_bytes == 1:
            dd_type = 'int8'

        if load_data:
            self.data = np.zeros((n_ints, n_ifs, n_chans_selected), dtype='float32')

            for ii in range(n_ints):
                """d = f.read(n_bytes * n_chans * n_ifs)
                """

                for jj in range(n_ifs):

                    f.seek(n_bytes * i0, 1) # 1 = from current location
                    #d = f.read(n_bytes * n_chans_selected)
                    #bytes_to_read = n_bytes * n_chans_selected

                    dd = np.fromfile(f, count=n_chans_selected, dtype=dd_type)

                    # Reverse array if frequency axis is flipped
                    if f_delt < 0:
                        dd = dd[::-1]

                    self.data[ii, jj] = dd

                    f.seek(n_bytes * (n_chans - i1), 1)  # Seek to start of next block
        else:
            print "Skipping data load..."
            self.data = np.array([0])

        ## Setup time axis
        t0 = self.header['tstart']
        t_delt = self.header['tsamp']
        self.timestamps = np.arange(0, n_ints) * t_delt / 24./60./60 + t0


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



    def write_file(self, filename_out):
        """ Write data to filterbank file.

        Args:
            filename_out (str): Name of output file
        """



        raise NotImplementedError('For now, writing files will be still implemented by filterbank.py')


        #rewrite header to be consistent with modified data
        self.header['fch1']   = self.freqs[0]
        self.header['foff']   = self.freqs[1] - self.freqs[0]
        self.header['nchans'] = self.freqs.shape[0]

        n_bytes  = self.__n_bytes

        with open(filename_out, "w") as fileh:
            fileh.write(generate_sigproc_header(self))
            j = self.data
            if n_bytes == 4:
                np.float32(j[:, ::-1].ravel()).tofile(fileh)
            elif n_bytes == 2:
                np.int16(j[:, ::-1].ravel()).tofile(fileh)
            elif n_bytes == 1:
                np.int8(j[:, ::-1].ravel()).tofile(fileh)


def open_file(filename,f_start=None, f_stop=None,t_start=None, t_stop=None):
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
        print filename
        raise IOError("No such directory: " + filename)

    filename = os.path.expandvars(os.path.expanduser(filename))
    # Get file extension to determine type
    ext = filename.split(".")[-1].strip().lower()


    if ext == 'h5':
        # Open HDF5 file
        return H5_reader(filename,f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)
    elif ext == 'fil':
        # Open FIL file
        return FIL_reader(filename,f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)
    else:
        # Fall back to regular Python `open` function
        return open(filename, *args, **kwargs)

