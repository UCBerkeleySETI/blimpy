import numpy as np

import logging
import psutil
import sys

logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    lformat = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    lformat = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'

logging.basicConfig(format=lformat, stream=stream, level=level_log)

# Convenient for memory calculations
GIGA = 1024 ** 3

# Threshold size for high resolution data
HIRES_THRESHOLD = 2**20


class Reader(object):
    """ Basic reader object """

    def __init__(self):

        self.t_begin = 0
        self.t_end   = 0

        # Calculate the max data array size from available memory
        self.available_memory = psutil.virtual_memory().available
        logger.debug("Reader __init__ available_memory={}".format(self.available_memory))
        if self.available_memory > GIGA:
            self.max_data_array_size = self.available_memory - GIGA
        else:
            self.max_data_array_size = self.available_memory
            logger.warning("Very low on memory, only {:.2f} MB available for use."
                           .format(float(self.available_memory) / 1e6))
        logger.debug("Reader __init__ max_data_array_size={}".format(self.max_data_array_size))

    def _setup_selection_range(self, f_start=None, f_stop=None, t_start=None, t_stop=None, init=False):
        """Making sure the selection if time and frequency are within the file limits.

        Args:
            init (bool): If call during __init__
        """


        def _dump_parms():
            wstr = "f_start={}, f_stop={}, t_start={}, t_stop={}, init={}" \
                  .format(f_start, f_stop, t_start, t_stop, init)
            logger.warning(wstr)


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
                _dump_parms()
            self.t_start = self.t_begin

        if t_stop <= self.t_end  and t_stop > self.t_begin:
            self.t_stop = int(t_stop)
        else:
            if init is False or t_stop:
                logger.warning('Setting t_stop = %f, since t_stop not given or not valid.'%self.t_end)
                _dump_parms()
            self.t_stop = self.t_end

        if f_start >= self.f_begin and f_start < self.f_end:
            self.f_start = f_start
        else:
            if init is False or f_start:
                logger.warning('Setting f_start = %f, since f_start not given or not valid.'%self.f_begin)
                _dump_parms()
            self.f_start = self.f_begin

        if f_stop <= self.f_end and f_stop > self.f_begin:
            self.f_stop = f_stop
        else:
            if init is False or f_stop:
                logger.warning('Setting f_stop = %f, since f_stop not given or not valid.'%self.f_end)
                _dump_parms()
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
            return 'float32'
        elif self._n_bytes  == 2:
            return 'uint16'
        elif self._n_bytes  == 1:
            return 'uint8'
        else:
            errmsg = 'Reader._setup_dtype detected invalid/unsupported header data: _n_bytes={}, nbits={}'.format(self._n_bytes, self.header['nbits'])
            logger.error(errmsg)
            raise ValueError(errmsg)

    def _calc_selection_size(self):
        """Calculate size of data of interest.
        """

        #Check to see how many integrations requested
        n_ints = self.t_stop - self.t_start
        #Check to see how many frequency channels requested
        n_chan = (self.f_stop - self.f_start) / abs(self.header['foff'])

        n_bytes  = self._n_bytes
        selection_size = int(n_ints*n_chan*n_bytes)

        return selection_size

    def _calc_selection_shape(self):
        """Calculate shape of data of interest.
        """

        #Check how many integrations requested
        n_ints = int(self.t_stop - self.t_start)
        #Check how many frequency channels requested
        n_chan = int(np.round((self.f_stop - self.f_start) / abs(self.header['foff'])))

        selection_shape = (n_ints,int(self.header['nifs']),n_chan)

        return selection_shape

    def _setup_chans(self):
        """Setup channel borders
        """

        if self.header['foff'] < 0:
            f0 = self.f_end
        else:
            f0 = self.f_begin

        i_start, i_stop = 0, self.n_channels_in_file
        if self.f_start:
            i_start = np.round((self.f_start - f0) / self.header['foff'])
        if self.f_stop:
            i_stop  = np.round((self.f_stop - f0)  / self.header['foff'])

        #calculate closest true index value
        chan_start_idx = int(i_start)
        chan_stop_idx  = int(i_stop)

        if chan_stop_idx < chan_start_idx:
            chan_stop_idx, chan_start_idx = chan_start_idx,chan_stop_idx

        self.chan_start_idx =  chan_start_idx
        self.chan_stop_idx = chan_stop_idx

    def _setup_freqs(self):
        """Updating frequency borders from channel values
        """

        if self.header['foff'] > 0:
            self.f_start = self.f_begin + self.chan_start_idx*abs(self.header['foff'])
            self.f_stop = self.f_begin + self.chan_stop_idx*abs(self.header['foff'])
        else:
            self.f_start = self.f_end - self.chan_stop_idx*abs(self.header['foff'])
            self.f_stop = self.f_end - self.chan_start_idx*abs(self.header['foff'])

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
        t0 = self.header['tstart']
        t_delt = self.header['tsamp']

        if update_header:
            timestamps = ii_start * t_delt / 24./60./60. + t0
        else:
            timestamps = np.arange(ii_start, ii_stop) * t_delt / 24./60./60. + t0

        return timestamps

    def populate_freqs(self):
        """
         Populate frequency axis
        """

        if self.header['foff'] < 0:
            f0 = self.f_end
        else:
            f0 = self.f_begin

        self._setup_chans()

        #create freq array
        i_vals = np.arange(self.chan_start_idx, self.chan_stop_idx)
        freqs = self.header['foff'] * i_vals + f0

        return freqs

    def adjust_n_coarse_chan(self, n_coarse_chan, nchans):
        r"""Don't let the calculated n_coarse_chan be < 1
            nor exceed the number of fine channels."""
        if n_coarse_chan < 1:
            errmsg1 = "blimpy:io:base_reader:adjust_n_coarse_chan: n_coarse_chan={}, nchans={}" \
                      .format(n_coarse_chan, nchans)
            logger.warning(errmsg1)
            errmsg2 = "blimpy:io:base_reader:adjust_n_coarse_chan: n_coarse_chan < 1. Replacing that with a value of 64 (SWAG)."
            logger.warning(errmsg2)
            return 64
        if n_coarse_chan > nchans: # exceeds the number of fine channels?
            errmsg1 = "blimpy:io:base_reader:adjust_n_coarse_chan: n_coarse_chan={}, nchans={}" \
                      .format(n_coarse_chan, nchans)
            logger.warning(errmsg1)
            errmsg2 = "blimpy:io:base_reader:adjust_n_coarse_chan: n_coarse_chan > nchans. Replacing that with the value of nchans (SWAG)."
            logger.warning(errmsg2)
            return nchans
        return n_coarse_chan

    def calc_n_coarse_chan(self, chan_bw=None):
        """ This makes an attempt to calculate the number of coarse channels in a given file.

            Note:
                This is unlikely to work on non-Breakthrough Listen data, as a-priori knowledge of
                the digitizer system is required.

        """

        nchans = int(self.header['nchans'])

        # Do we have a file with enough channels that it has coarse channelization?
        if chan_bw is not None:
            bandwidth = abs(self.f_stop - self.f_start)
            n_coarse_chan = bandwidth / chan_bw
            return self.adjust_n_coarse_chan(n_coarse_chan, nchans)

        # High resolution data?
        if nchans >= HIRES_THRESHOLD:
            # Does the common FFT length of 2^20 divide through without a remainder?
            # This should work for most GBT and all Parkes hires data
            if nchans % HIRES_THRESHOLD == 0:
                n_coarse_chan = nchans / float(HIRES_THRESHOLD)
                return self.adjust_n_coarse_chan(n_coarse_chan, nchans)
            # Early GBT data has non-2^N FFT length, check if it is GBT data
            elif self.header['telescope_id'] == 6:
                coarse_chan_bw = 2.9296875
                bandwidth = abs(self.f_stop - self.f_start)
                n_coarse_chan = bandwidth / coarse_chan_bw
                return self.adjust_n_coarse_chan(n_coarse_chan, nchans)
            else:
                errmsg1 = "blimpy:io:base_reader:calc_n_coarse_chan: hires nchans not divisible by 2^20 and not GBT"
                logger.warning(errmsg1)
                errmsg2 = "Setting a value of 64 (SWAG). In turbo_seti, you can specify n_course_chan explicitly."
                logger.info(errmsg2)
                return 64
                
        # Not high resolution data.  GBT?
        if self.header['telescope_id'] == 6:
            #For GBT non-hires data
            coarse_chan_bw = 2.9296875
            bandwidth = abs(self.f_stop - self.f_start)
            n_coarse_chan = bandwidth / coarse_chan_bw
            return self.adjust_n_coarse_chan(n_coarse_chan, nchans)
            
        # Not high resolution data and Not GBT
        else: 
            errmsg1 = "blimpy:io:base_reader:calc_n_coarse_chan: not hires and not GBT. Setting a value of 1"
            logger.warning(errmsg1)
            errmsg2 = "Setting a value of 64 (SWAG). In turbo_seti, you can specify n_course_chan explicitly."
            logger.info(errmsg2)
            return 64

    def calc_n_blobs(self, blob_dim):
        """ Given the blob dimensions, calculate how many fit in the data selection.
        """

        n_blobs = int(np.ceil(1.0 * np.prod(self.selection_shape) / np.prod(blob_dim)))

        return n_blobs

    def warn_memory(self, name, size):
        """ Warn that <name> is larger than our limit for loading things into memory.
        """
        logger.warning(f"{name} size is {size / GIGA:.2f} GB, which exceeds the memory usage limit of {self.max_data_array_size / GIGA} GB. Keeping data on disk.")
    
    def isheavy(self):
        """ Check if the current selection is too large.
        """

        selection_size_bytes = self._calc_selection_size()

        if selection_size_bytes > self.max_data_array_size:
            return True
        else:
            return False
