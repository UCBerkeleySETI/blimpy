#!/usr/bin/env python
"""
# waterfall.py

Python class and command line utility for reading and plotting waterfall files.

This provides a class, Waterfall(), which can be used to read a blimpy file (.fil or .h5):


    fil = Waterfall("test_psr.fil")
    print(fil.header)
    print(fil.data.shape)
    print(fil.freqs)

    plt.figure()
    fil.plot_spectrum(t=0)
    plt.show()

"""

import sys
import os
import numpy as np
import six

SHOWING_BACKEND = False

from blimpy.io import file_wrapper as fw
from .plotting.config import plt, get_mpl_backend, set_mpl_backend, print_plotting_backend, ok_to_show
from .plotting import plot_all, plot_kurtosis, plot_spectrum_min_max, plot_spectrum, plot_time_series, plot_waterfall
MPL_BACKEND = get_mpl_backend()
if SHOWING_BACKEND:
    print(f"after importing plot_all etc: {MPL_BACKEND}")

from astropy.time import Time
from astropy import units as u

# Logging set up
import logging
logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    fmt = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    fmt = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'

logging.basicConfig(format=fmt, stream=stream, level=level_log)

MAX_BLOB_MB = 1024


###
# Main blimpy class
###

class Waterfall():
    """ Class for loading and writing blimpy data (.fil, .h5) """


    def __repr__(self):
        return "Waterfall data: %s" % self.filename


    def _init_alternate(self, header_dict, data, filename=None):

        # Validate parameters.
        assert filename is None
        assert isinstance(header_dict, dict)
        assert "nchans" in header_dict
        assert "fch1" in header_dict
        assert data.ndim == 3
        assert data.shape[1] == 1

        # Set dummy-file Waterfall object properties.
        self.filename = None
        self.ext = self.filename
        self.file_shape = self.filename
        self.file_size_bytes = 0

        # The Waterfall "container":
        class Container():
            n_beams_in_file = 1
            n_pols_in_file = 1
            _d_type = np.float32
            t_begin = 0
            freq_axis = 2
            time_axis = 0
            beam_axis = 1
        self.container = Container()
        self.container.n_channels_in_file  = header_dict["nchans"]
        self.container._n_bytes = int(header_dict["nbits"] / 8)  # number of bytes per digit.
        if header_dict['foff'] < 0:
            self.container.f_end  = header_dict['fch1']
            self.container.f_begin  = self.container.f_end + self.container.n_channels_in_file * header_dict['foff']
        else:
            self.container.f_begin  = header_dict['fch1']
            self.container.f_end  = self.container.f_begin + self.container.n_channels_in_file * header_dict['foff']
        self.container.f_start = self.container.f_begin
        self.container.f_stop = self.container.f_end
        self.container.t_end = data.shape[0]

        # Set Waterfall object properties.
        self.header = header_dict
        self.file_header = header_dict
        self.n_ints_in_file = data.shape[0]
        self.selection_shape = data.shape
        self.n_channels_in_file = header_dict["nchans"]
        self.freq_axis = 2
        self.time_axis = 0
        self.beam_axis = 1
        self.stokes_axis = 4
        self.logger = logger

        # Attach data matrix.
        self.data = data

        # Attach plotting methods.
        self.plot_spectrum         = six.create_bound_method(plot_spectrum, self)
        self.plot_waterfall        = six.create_bound_method(plot_waterfall, self)
        self.plot_kurtosis         = six.create_bound_method(plot_kurtosis, self)
        self.plot_time_series      = six.create_bound_method(plot_time_series, self)
        self.plot_all              = six.create_bound_method(plot_all, self)
        self.plot_spectrum_min_max = six.create_bound_method(plot_spectrum_min_max, self)


    def __init__(self, filename=None, f_start=None, f_stop=None, t_start=None, t_stop=None,
                 load_data=True, max_load=None, header_dict=None, data_array=None):
        """ Class for loading and plotting blimpy data.

        This class parses the blimpy file and stores the header and data
        as objects:
            fb = Waterfall('filename_here.fil')
            fb.header        # blimpy header, as a dictionary
            fb.data          # blimpy data, as a numpy array

        Args:
            filename (str): filename of blimpy file.  REQUIRED.
            f_start (float): start frequency in MHz
            f_stop (float): stop frequency in MHz
            t_start (int): start integration ID
            t_stop (int): stop integration ID
            load_data (bool): load data. If set to False, only header will be read.
            max_load (float): maximum data to load in GB.
            header_dict (dict): *NOT CURRENTLY SUPPORTED*
            data_array (np.array): *NOT CURRENTLY SUPPORTED*
        """

        if (header_dict is not None) or (data_array is not None):
            self._init_alternate(header_dict, data_array, filename=filename)
            return

        if filename is None:
            raise ValueError("Currently, a value for filename must be supplied.")

        self.filename = filename
        self.ext = os.path.splitext(filename)[-1].lower()
        self.container = fw.open_file(filename, f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop,
                                      load_data=load_data, max_load=max_load)
        self.file_header = self.container.header
        self.header = self.file_header
        self.n_ints_in_file = self.container.n_ints_in_file
        self.file_shape = self.container.file_shape
        self.file_size_bytes = self.container.file_size_bytes
        self.selection_shape = self.container.selection_shape
        self.n_channels_in_file = self.container.n_channels_in_file

        # These values will be modified once code for multi_beam and multi_stokes observations are possible.
        self.freq_axis = 2
        self.time_axis = 0
        self.beam_axis = 1  # Place holder # Polarisation?
        self.stokes_axis = 4  # Place holder

        self.logger = logger

        self.__load_data()

        # Attach methods
        self.plot_spectrum         = six.create_bound_method(plot_spectrum, self)
        self.plot_waterfall        = six.create_bound_method(plot_waterfall, self)
        self.plot_kurtosis         = six.create_bound_method(plot_kurtosis, self)
        self.plot_time_series      = six.create_bound_method(plot_time_series, self)
        self.plot_all              = six.create_bound_method(plot_all, self)
        self.plot_spectrum_min_max = six.create_bound_method(plot_spectrum_min_max, self)


    def __load_data(self):
        """ Helper for loading data from a container. Should not be called manually. """

        self.data = self.container.data

    def get_freqs(self):
        """
        Get the frequency array for this Waterfall object.

        Returns
        -------
        numpy array
            Values for all of the fine frequency channels.

        """
        return np.arange(0, self.header['nchans'], 1, dtype=float) \
                          * self.header['foff'] + self.header['fch1']

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None):
        """ Reads data selection if small enough.

        Args:
            f_start (float): Start frequency in MHz
            f_stop (float): Stop frequency in MHz
            t_start (int): Integer time index to start at
            t_stop (int): Integer time index to stop at

        Data is loaded into self.data (nothing is returned)
        """

        self.container.read_data(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)

        self.__load_data()

    def _update_header(self):
        """ Updates the header information from the original file to the selection. """

        #Updating frequency of first channel from selection
        if self.header['foff'] < 0:
            self.header['fch1'] = self.container.f_stop
        else:
            self.header['fch1'] = self.container.f_start

        #Updating number of fine channels.
        self.header['nchans'] = self.container.selection_shape[self.freq_axis]

        #Updating time stamp for first time bin from selection
        self.header['tstart'] = self.container.populate_timestamps(update_header=True)

    def info(self):
        """ Print header information and other derived information. """

        print("\n--- File Info ---")

        for key, val in self.file_header.items():
            if key == 'src_raj':
                val = val.to_string(unit=u.hour, sep=':')
            if key == 'src_dej':
                val = val.to_string(unit=u.deg, sep=':')
            if key in ('foff', 'fch1'):
                val *= u.MHz
            if key == 'tstart':
                print("%16s : %32s" % ("tstart (ISOT)", Time(val, format='mjd').isot))
                key = "tstart (MJD)"
            print("%16s : %32s" % (key, val))

        print("\n%16s : %32s" % ("Num ints in file", self.n_ints_in_file))
        print("%16s : %32s" % ("File shape", self.file_shape))
        print("--- Selection Info ---")
        print("%16s : %32s" % ("Data selection shape", self.selection_shape))
        if self.header['foff'] < 0: # descending frequency values
            minfreq = self.container.f_start - self.header['foff']
            maxfreq = self.container.f_stop
        else: # ascending frequency values
            minfreq = self.container.f_start
            maxfreq = self.container.f_stop - self.header['foff']
        print("%16s : %32s" % ("Minimum freq (MHz)", minfreq))
        print("%16s : %32s" % ("Maximum freq (MHz)", maxfreq))


    def _get_blob_dimensions(self, chunk_dim):
        """ Computes the blob dimensions, trying to read around 1024 MiB at a time.
            This is assuming a chunk is about 1 MiB.

            Notes:
                A 'blob' is the max size that will be read into memory at once.
                A 'chunk' is a HDF5 concept to do with efficient read access, see
                https://portal.hdfgroup.org/display/HDF5/Chunking+in+HDF5

            Args:
                chunk_dim (array of ints): Shape of chunk, e.g. (1024, 1, 768)

            Returns blob dimensions (shape = array of ints).
        """

        #Taking the size into consideration, but avoiding having multiple blobs within a single time bin.
        freq_axis_size = self.selection_shape[self.freq_axis]
        if self.selection_shape[self.freq_axis] > chunk_dim[self.freq_axis] * MAX_BLOB_MB:
            time_axis_size = 1
        else:
            time_axis_size = np.min([chunk_dim[self.time_axis] * MAX_BLOB_MB * chunk_dim[self.freq_axis] / freq_axis_size, self.selection_shape[self.time_axis]])

        blob_dim = (int(time_axis_size), 1, freq_axis_size)

        return blob_dim

    def _get_chunk_dimensions(self):
        """ Sets the chunking dimensions depending on the file type.

            Notes: A 'chunk' is a HDF5 concept to do with efficient read access, see
            https://portal.hdfgroup.org/display/HDF5/Chunking+in+HDF5

            Returns chunk dimensions (shape), e.g. (2048, 1, 512)
        """

        #Usually '.0000.' is in self.filename
        if np.abs(self.header['foff']) < 1e-5:
            logger.info('Detecting high frequency resolution data.')
            chunk_dim = (1,1,1048576) #1048576 is the number of channels in a coarse channel.
            return chunk_dim

        #Usually '.0001.' is in self.filename
        if np.abs(self.header['tsamp']) < 1e-3:
            logger.info('Detecting high time resolution data.')
            chunk_dim = (2048,1,512) #512 is the total number of channels per single band (ie. blc00)
            return chunk_dim

        #Usually '.0002.' is in self.filename
        if np.abs(self.header['foff']) < 1e-2 and np.abs(self.header['foff'])  >= 1e-5:
            logger.info('Detecting intermediate frequency and time resolution data.')
            chunk_dim = (10,1,65536)  #65536 is the total number of channels per single band (ie. blc00)
            return chunk_dim

        logger.warning('File format not known. Will use minimum chunking. NOT OPTIMAL.')
        chunk_dim = (1,1,512)
        return chunk_dim

    def calc_n_coarse_chan(self, chan_bw=None):
        """ This makes an attempt to calculate the number of coarse channels in a given freq selection.

            Note:
                This is unlikely to work on non-Breakthrough Listen data, as a-priori knowledge of
                the digitizer system is required.

            Returns n_coarse_chan (int), number of coarse channels
        """

        n_coarse_chan = self.container.calc_n_coarse_chan(chan_bw)

        return n_coarse_chan

    def grab_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None, if_id=0):
        """ Extract a portion of data by frequency range.

        Args:
            f_start (float): start frequency in MHz
            f_stop (float): stop frequency in MHz
            if_id (int): IF input identification (req. when multiple IFs in file)

        Returns:
            (freqs, data) (np.arrays): frequency axis in MHz and data subset
        """

        if self.container.isheavy():
            raise Exception("Waterfall.grab_data: Large data array was not loaded!  Try instantiating Waterfall with max_load.")

        try:
            self.freqs = self.container.populate_freqs()
        except Exception as ex:
            raise Exception("Waterfall.grab_data: Cannot allocate frequency array") from ex

        try:
            self.timestamps = self.container.populate_timestamps()
        except Exception as ex:
            raise Exception("Waterfall.grab_data: Cannot allocate timestamps array") from ex

        if f_start is None:
            f_start = self.freqs[0]
        if f_stop is None:
            f_stop = self.freqs[-1]

        try:
            i0 = np.argmin(np.abs(self.freqs - f_start))
            i1 = np.argmin(np.abs(self.freqs - f_stop))

            if i0 < i1:
                plot_f    = self.freqs[i0:i1 + 1]
                plot_data = np.squeeze(self.data[t_start:t_stop, if_id, i0:i1 + 1])
            else:
                plot_f    = self.freqs[i1:i0 + 1]
                plot_data = np.squeeze(self.data[t_start:t_stop, if_id, i1:i0 + 1])
        except:
            raise Exception("Waterfall.grab_data: Too much data requested")

        return plot_f, plot_data

    def write_to_fil(self, filename_out, *args, **kwargs):
        from blimpy.io import write_to_fil
        write_to_fil(self, filename_out)

    def write_to_hdf5(self, filename_out, *args, **kwargs):
        from blimpy.io import write_to_hdf5
        write_to_hdf5(self, filename_out, *args, **kwargs)

    def blank_dc(self, n_coarse_chan):
        """ Blank DC bins in coarse channels.

        Removes the DC spike in centre of coarse channel bins.

        Note: currently only works if entire file is read
        """

        if n_coarse_chan < 1:
            logger.warning('Coarse channel number < 1, unable to blank DC bin.')
            return

        if not n_coarse_chan % int(n_coarse_chan) == 0:
            logger.warning('Selection does not contain an integer number of coarse channels, unable to blank DC bin.')
            return

        n_coarse_chan = int(n_coarse_chan)

        n_chan = self.data.shape[-1]
        n_chan_per_coarse = int(n_chan / n_coarse_chan) # ratio of fine channels to coarse channels

        mid_chan = int(n_chan_per_coarse / 2)

        def parse(data, n_coarse_chan, n_chan_per_coarse, mid_chan):
            for ii in range(n_coarse_chan):
                ss = ii*n_chan_per_coarse
                w_slice = data[..., ss+mid_chan+5:ss+mid_chan+10]
                chtest = w_slice.shape[-1]
                # If we are nearing the end of the fine channel frequency array,
                # do not do anymore.  See issue #212.
                if chtest < 5:
                    break
                data[..., ss+mid_chan] = np.median(w_slice)
        
        parse(self.data, n_coarse_chan, n_chan_per_coarse, mid_chan)

    def calibrate_band_pass_N1(self):
        """ One way to calibrate the band pass is to take the median value
            for every frequency fine channel, and divide by it.

            sets data = data / band_pass
        """

        band_pass = np.median(self.data.squeeze(),axis=0)
        self.data = self.data/band_pass


def cmd_tool(args=None):
    """ Command line tool for plotting and viewing info on blimpy files """


    from argparse import ArgumentParser

    set_mpl_backend(MPL_BACKEND)
    if SHOWING_BACKEND:
        print_plotting_backend("entered cmd_tool")

    parser = ArgumentParser(description="Command line utility for reading and plotting blimpy files.")

    parser.add_argument('filename', type=str,
                        help='Name of file to read')
    parser.add_argument('-p', action='store', dest='what_to_plot', type=str,
                        choices=['w', 's', 't', 'k', 'mm', 'a', 'ank'],
                        help='Show: "w" waterfall (freq vs. time) plot; "s" integrated spectrum plot; \
                        "t" for time series; "mm" for spectrum including min max; "k" for kurtosis; \
                        "a" for all available plots and information; and "ank" for all but kurtosis.')
    parser.add_argument('-b', action='store', default=None, dest='f_start', type=float,
                        help='Start frequency (begin), in MHz')
    parser.add_argument('-e', action='store', default=None, dest='f_stop', type=float,
                        help='Stop frequency (end), in MHz')
    parser.add_argument('-B', action='store', default=None, dest='t_start', type=int,
                        help='Start integration (begin, inclusive) ID ')
    parser.add_argument('-E', action='store', default=None, dest='t_stop', type=int,
                        help='Stop integration (end, exclusive) ID')
    parser.add_argument('-i', action='store_true', default=False, dest='info_only',
                        help='Show info only')
    parser.add_argument('-a', action='store_true', default=False, dest='average',
                       help='average along time axis (plot spectrum only)')
    parser.add_argument('-s', action='store', default='', dest='plt_filename', type=str,
                       help='save plot graphic to file (give filename as argument)')
    parser.add_argument('-S', action='store_true', default=False, dest='save_only',
                       help='Turn off plotting of data and only save to file.')
    parser.add_argument('-D', action='store_false', default=False, dest='blank_dc',
                       help='Use to blank DC bin.')
    parser.add_argument('-H', action='store_true', default=False, dest='to_hdf5',
                       help='Write file to hdf5 format.')
    parser.add_argument('-F', action='store_true', default=False, dest='to_fil',
                       help='Write file to .fil format.')
    parser.add_argument('-o', action='store', default=None, dest='filename_out', type=str,
                        help='Filename output (if not provided, the name will be the same but with appropriate extension).')
    parser.add_argument('-l', action='store', default=None, dest='max_load', type=float,
                        help='Maximum data limit to load. Default:1GB')

    if args is None:
        args = sys.argv[1:]

    parse_args = parser.parse_args(args)

    # Open blimpy data
    filename = parse_args.filename
    load_data = not parse_args.info_only
    filename_out = parse_args.filename_out

    fil = Waterfall(filename,
                    f_start=parse_args.f_start,
                    f_stop=parse_args.f_stop,
                    t_start=parse_args.t_start,
                    t_stop=parse_args.t_stop,
                    load_data=load_data,
                    max_load=parse_args.max_load)
    fil.info()
    if parse_args.info_only:
        return

    # Blank DC.

    if parse_args.blank_dc:
        logger.info("Blanking DC bin")
        n_coarse_chan = fil.calc_n_coarse_chan()
        fil.blank_dc(n_coarse_chan)

    # Plotting

    if parse_args.what_to_plot is not None:

        if SHOWING_BACKEND:
            print_plotting_backend("before plotting")

        if parse_args.what_to_plot == "s":
            plt.figure("Spectrum", figsize=(8, 6))
            fil.plot_spectrum(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all')
        elif parse_args.what_to_plot == "mm":
            plt.figure("min max", figsize=(8, 6))
            fil.plot_spectrum_min_max(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all')
        elif parse_args.what_to_plot == "k":
            plt.figure("kurtosis", figsize=(8, 6))
            fil.plot_kurtosis(f_start=parse_args.f_start, f_stop=parse_args.f_stop)
        elif parse_args.what_to_plot == "t":
            plt.figure("Time Series", figsize=(8, 6))
            fil.plot_time_series(f_start=parse_args.f_start, f_stop=parse_args.f_stop,orientation='h')
        elif parse_args.what_to_plot == "a":
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all')
        elif parse_args.what_to_plot == "ank":
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all', kurtosis=False)
        else: # parse_args.what_to_plot = "w"
            plt.figure("waterfall", figsize=(8, 6))
            fil.plot_waterfall(f_start=parse_args.f_start, f_stop=parse_args.f_stop)

        if parse_args.plt_filename != '':
            plt.savefig(parse_args.plt_filename)

        if SHOWING_BACKEND:
            print_plotting_backend("before plt.show")
        if not parse_args.save_only:
            if ok_to_show():
                plt.show()
            else:
                logger.warning("No $DISPLAY available.")

    # Save in specified Filterbank format

    fileroot = os.path.splitext(filename)[0]

    if parse_args.to_hdf5 and parse_args.to_fil:
        raise ValueError('Either provide to_hdf5 or to_fil, but not both.')

    if parse_args.to_hdf5:
        if not filename_out:
            filename_out = fileroot + '.h5'

        logger.info(f"Writing FBH5 file : {filename_out}")
        fil.write_to_hdf5(filename_out)
        logger.info('File written.')

    if parse_args.to_fil:
        if not filename_out:
            filename_out = fileroot + '.fil'

        logger.info(f"Writing SIGPROC Filterbank file : {filename_out}")
        fil.write_to_fil(filename_out)
        logger.info('File written.')


if __name__ == "__main__":
    cmd_tool()
