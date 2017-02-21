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
MAX_DATA_ARRAY_SIZE = 1024 * 1024 * 1024. * 8    # Max size of data array to load into memory (in bytes)
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
            self.n_ints_in_file  = self.h5["data"].shape[0] #
            self.n_channels_in_file  = self.h5["data"].shape[2] #
            self.n_beams_in_file = 1 #Placeholder for future development.
            self.n_pols_in_file = 0 #Placeholder for future development.
            self.data_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)
            self.__setup_time_axis()

            if self.file_size_bytes > MAX_DATA_ARRAY_SIZE:
                self.heavy = True
                if f_start or f_stop or t_start or t_stop:
                    selection_size = self.__calc_selection_size(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)
                    if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
                        logger.warning("Selection size of %f MB, exceeding our size limit %f MB. Data not loaded, please try another (t,v) selection."%(self.selection_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                    else:
                        self.read_data(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)
                else:
                    logger.warning("The file is of size %f MB, exceeding our size limit %f MB. Data not loaded."%(self.file_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            else:
                self.heavy = False
                self.read_data(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)

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

    def __calc_selection_size(f_start=None, f_stop=None, t_start=None, t_stop=None):
        '''Calculate size of data of interest.
        '''

        #Check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if t_start:
            ii_start = t_start
        if t_stop:
            ii_stop = t_stop
        n_ints = ii_stop - ii_start

        #Check to see how many frequency channels requested
        jj_start, jj_stop = 0, self.n_channels_in_file
        if f_start:
            jj_start = f_start
        if f_stop:
            jj_stop = f_stop
        n_chan = (jj_stop - jj_start) / abs(self.header['foff'])

        selection_size = n_ints*n_chan*32/(8.)

        return selection_size

    def read_data(self, f_start=None, f_stop=None, t_start=None, t_stop=None, load_data=True):

        self.data = self.h5["data"][:]
        self.__setup_freqs()

    def __setup_freqs(self, f_start=None, f_stop=None):
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
            self.__get_n_ints_in_file()
            self.__setup_time_axis()
            self.n_channels_in_file  = self.header['nchans']
            self.n_beams_in_file = 1 #Placeholder for future development.
            self.n_pols_in_file = 0 #Placeholder for future development.
            self.data_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)

            # set start of data, at real length of header  (future development.)
#            self.datastart=self.hdrraw.find('HEADER_END')+len('HEADER_END')+self.startsample*self.channels

            if self.file_size_bytes > MAX_DATA_ARRAY_SIZE:
                self.heavy = True
                if f_start or f_stop or t_start or t_stop:
                    selection_size = self.__calc_selection_size(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)
                    if selection_size_bytes > MAX_DATA_ARRAY_SIZE:
                        logger.warning("Selection size of %f MB, exceeding our size limit %f MB. Data not loaded, please try another (t,v) selection."%(self.selection_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
                    else:
                        self.read_data(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)
                else:
                    logger.warning("The file is of size %f MB, exceeding our size limit %f MB. Data not loaded."%(self.file_size_bytes/(1024.**2), MAX_DATA_ARRAY_SIZE/(1024.**2)))
            else:
                self.heavy = False
                self.read_data(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)

        else:
            raise IOError("Need a file to open, please give me one!")

    def __calc_selection_size(f_start=None, f_stop=None, t_start=None, t_stop=None):
        '''Calculate size of data of interest.
        '''

        #Check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
        if t_start:
            ii_start = t_start
        if t_stop:
            ii_stop = t_stop
        n_ints = ii_stop - ii_start

        #Check to see how many frequency channels requested
        jj_start, jj_stop = 0, self.n_channels_in_file
        if f_start:
            jj_start = f_start
        if f_stop:
            jj_stop = f_stop
        n_chan = (jj_stop - jj_start) / abs(self.header['foff'])

        selection_size = n_ints*n_chan*32/(8.)

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

        n_bytes  = self.header['nbits'] / 8
        n_chans = self.header['nchans']
        n_ifs   = self.header['nifs']

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

    def __setup_freqs(self, f_start=None, f_stop=None):
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

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True):

        ## Setup frequency axis
        f0 = self.header['fch1']
        f_delt = self.header['foff']

        # keep this seperate!
        # file_freq_mapping =  np.arange(0, self.header['nchans'], 1, dtype='float64') * f_delt + f0

        #convert input frequencies into what their corresponding index would be

        i_start, i_stop, chan_start_idx, chan_stop_idx = self.__setup_freqs(f_start, f_stop)

        n_bytes  = self.header['nbits'] / 8
        n_chans = self.header['nchans']
        n_chans_selected = self.freqs.shape[0]
        n_ifs   = self.header['nifs']

        # Load binary data
        f = open(self.filename, 'rb')
        f.seek(self.idx_data)

        # now check to see how many integrations requested
        ii_start, ii_stop = 0, self.n_ints_in_file
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

        #Set up the data type (taken out of loop for speed)
        if n_bytes == 4:
            dd_type = 'float32'
        elif n_bytes == 2:
            dd_type = 'int16'
        elif n_bytes == 1:
            dd_type = 'int8'

        if load_data:

            #EE added n_bytes here, lets thing again about this...
            if n_ints * n_ifs * n_chans_selected*n_bytes > MAX_DATA_ARRAY_SIZE:
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


def open_file(filename, *args, **kwargs):
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
        return H5_reader(filename, *args, **kwargs)
    elif ext == 'fil':
        # Open FIL file
        return FIL_reader(filename, *args, **kwargs)
    else:
        # Fall back to regular Python `open` function
        return open(filename, *args, **kwargs)


class  nonFITS:
    """ This class is where the filterbank data is loaded, as well as the header info.
        It creates other atributes related to the search (load_drift_indexes).
        Similar to FITS, but in this case to load fil not fits.

    """
    def __init__(self, filename=None, size_limit = 1024.0):
        if filename and os.path.isfile(filename):
            self.filename = filename
            self.filestat = os.stat(filename)
            filesize = self.filestat.st_size/(1024.0**2)
            if filesize > size_limit:
                logger.error("The file is of size %f MB, exceeding our size limit %f MB. Aborting..."%(filesize, size_limit))
                return None
            try:
                fil_file2=fr2.DataReader(filename)  # Will be replaced by danny's filterbank...
                fil_file=fr.Filterbank(filename)
                header = self.make_fits_header(fil_file2.headerinfo)
#EE_fil2                header = self.make_fits_header(fil_file.header)
            except:
                logger.error("Error encountered when trying to open FITS file %s"%filename)
                self.status = False
                return None
            self.fftlen = header['NAXIS1']
            self.tsteps_valid = header['NAXIS2']
            self.tsteps = int(math.pow(2, math.ceil(np.log2(math.floor(self.tsteps_valid)))))   ## what is this for??
            self.obs_length = self.tsteps_valid * header['DELTAT']
            self.shoulder_size = 0
            self.tdwidth = self.fftlen + self.shoulder_size*self.tsteps  ##EE why is this multiplied by 8? This gives two regions, each of 4*steps, around spectra[i]
            self.drift_rate_resolution = (1e6 * header['DELTAF']) / self.obs_length
            self.nom_max_drift = self.drift_rate_resolution * self.tsteps_valid  ##EE Do I need tsteps here?

            ##EE: debug. Skyping barycenter for now. Would need to debug it first.
            if logger.getEffectiveLevel() > 1000:  ##EE: logging.getLevelName(10)='DEBUG'

                self.header = barycenter.correct(header, self.obs_length)
                logger.info('barycenter done for fits file %s! baryv: %f'%(filename, self.header['baryv']))
            else:
                self.header = header
                self.header['baryv'] = 0.0
                self.header['barya'] = 0.0

            # some default values
            self.original_vals= {'tsteps_valid': self.tsteps_valid, 'tsteps': self.tsteps,
                                 'tdwidth': self.tdwidth, 'fftlen':self.fftlen}
            self.compressed_t = False
            self.compressed_f = False
            self.status = True

    @staticmethod
    def make_fits_header(header,LOFAR=False):
        '''Takes .fil header into fits header format '''

        base_header = {}
        base_header['SIMPLE'] = True
        base_header['NAXIS'] = 2

        base_header['DOPPLER'] = 0.0
        base_header['SNR'] = 0.0
        base_header['EXTEND'] = True
        base_header['XTENSION'] = 'IMAGE   '
        base_header['PCOUNT'] = 1
        base_header['GCOUNT'] = 1

        if '32' in header['Number of bits per sample']:
            base_header['BITPIX'] = -32
        else:
            raise ValueError('Check nbits per sample. Not equeal 32')

        base_header['NAXIS1'] = int(header['Number of channels'])  #nchans
#EE_fil2        base_header['NAXIS1'] = int(header['nchans'])  #nchans
        base_header['NAXIS2'] = int(header['Number of samples'])
#EE_fil2        base_header['NAXIS2'] = int(header[''])
        base_header['DELTAT'] = float(header['Sample time (us)'])/1e6
        base_header['MJD'] = float(header['Time stamp of first sample (MJD)'])
        base_header['TOFFSET'] = float(header['Sample time (us)'])/1e6
#EE_fil2        base_header['DELTAT'] = float(header['tsamp'])
#EE_fil2        base_header['MJD'] = float(header['tstart'])
#EE_fil2        base_header['TOFFSET'] = float(header['tsamp)'])
        base_header['DELTAF'] =  np.abs(float(header['Channel bandwidth      (MHz)']))
#EE_fil2        base_header['DELTAF'] =  np.abs(float(header['foff']))
        base_header['SOURCE'] = header['Source Name'].replace('\xc2\xa0','_').replace(' ','')  #Removing white spaces and bad formats
#EE_fil2        base_header['SOURCE'] = header['source_name'].replace('\xc2\xa0','_').replace(' ','')   #Removing white spaces and bad formats
        base_header['FCNTR'] = float(header['Frequency of channel 1 (MHz)']) - base_header['DELTAF']*base_header['NAXIS1']/2
#EE_fil2        base_header['FCNTR'] = float(header['fch1']) - base_header['DELTAF']*base_header['NAXIS1']/2
        base_header['DEC'] = float(header['Source DEC (J2000)'])
        base_header['RA'] = float(header['Source RA (J2000)'])
#EE_fil2        base_header['DEC'] = float(header['Source DEC (J2000)'])
#EE_fil2        base_header['RA'] = float(header['Source RA (J2000)'])

        return base_header

    def load_data(self, max_search_rate=None, bw_compress_width=None, logwriter=None):
        ''' Read all the data from file.
        '''

#EE_fil        fil_file = fr2.DataReader(self.filename)
        fil_file = fr.Filterbank(filename)
#EE_fil        spec = fil_file.read_all()
        spec = np.squeeze(fil_file.data)
        spectra = np.array(spec, dtype=np.float64)

        if spectra.shape != (self.tsteps_valid, self.fftlen):
            raise ValueError('Something is wrong with array size.')

        drift_indexes = self.load_drift_indexes()

        return spectra, drift_indexes

    def load_drift_indexes(self):
        ''' The drift indexes are read from an stored file so that no need to recalculate. This speed things up.
        '''
        n = int(np.log2(self.tsteps))
        if n > 9:
            di_array = np.genfromtxt(resource_filename('dedoppler_bones', '../drift_indexes/drift_indexes_array_%d.txt'%n), delimiter=' ', dtype=int)
        else:
            di_array = np.genfromtxt(resource_filename('dedoppler_bones', '../drift_indexes/drift_indexes_array_%d.txt'%n), delimiter='\t', dtype=int)

        ts2 = self.tsteps/2
        drift_indexes = di_array[self.tsteps_valid - 1 - ts2, 0:self.tsteps_valid]
        return drift_indexes

    def get_info(self):
        return ""
