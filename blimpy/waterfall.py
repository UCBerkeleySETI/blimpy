#!/usr/bin/env python
"""
# waterfall.py

Python class and command line utility for reading and plotting waterfall files.

This provides a class, Waterfall(), which can be used to read a blimpy file (.fil or .h5):

    ````
    fil = Waterfall('test_psr.fil')
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

from .filterbank import Filterbank
from . import file_wrapper as fw
from .sigproc import *

try:
    import h5py
    HAS_HDF5 = True
except ImportError:
    HAS_HDF5 = False

try:
    HAS_BITSHUFFLE = True
    import bitshuffle.h5
except ImportError:
    HAS_BITSHUFFLE = False
    pass

#import pdb #pdb.set_trace()

# Check if $DISPLAY is set (for handling plotting on remote machines with no X-forwarding)
if 'DISPLAY' in os.environ.keys():
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
MAX_BLOB_MB         = 1024                    # Max size of blob in MB




###
# Main blimpy class
###

class Waterfall(Filterbank):
    """ Class for loading and writing blimpy data (.fil, .h5) """

    def __init__(self, filename=None, f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True,header_dict=None, data_array=None):
        """ Class for loading and plotting blimpy data.

        This class parses the blimpy file and stores the header and data
        as objects:
            fb = Waterfall('filename_here.fil')
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
        """

##EE        super(Waterfall, self).__init__()

        if filename:
            self.filename = filename
            self.ext = filename.split(".")[-1].strip().lower()  #File extension
            self.container = fw.open_file(filename, f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop,load_data=load_data)
            self.header = self.container.header
            self.n_ints_in_file = self.container.n_ints_in_file
            self.container.setup_time_axis()
            self.file_shape = self.container.file_shape
            self.file_size_bytes = self.container.file_size_bytes
            self.selection_shape = self.container.selection_shape
            self.n_channels_in_file = self.container.n_channels_in_file

            # These values will be modified once code for multi_beam and multi_stokes observations are possible.
            self.freq_axis = 2
            self.time_axis = 0
            self.beam_axis = 1  # Place holder
            self.stokes_axis = 4  # Place holder

            self.__load_data()

        elif header_dict is not None and data_array is not None:
            self.gen_from_header(header_dict, data_array)
        else:
            pass


    def __load_data(self):
        """
        """

        self.data = self.container.data
        self.freqs = self.container.freqs
        self.timestamps = self.container.timestamps

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None):
        """ Reads data selection if small enough.
        """

        self.container.read_data(f_start=f_start, f_stop=f_stop,t_start=t_start, t_stop=t_stop)

        self.__load_data()


    def write_to_fil(self, filename_out, *args, **kwargs):
        """ Write data to .fil file.
            It check the file size then decides how to write the file.

        Args:
            filename_out (str): Name of output file
        """

        if self.container.isheavy():
            self.__write_to_fil_heavy(filename_out)
        else:
            self.__write_to_fil_light(filename_out)

    def __write_to_fil_heavy(self, filename_out, *args, **kwargs):
        """ Write data to .fil file.

        Args:
            filename_out (str): Name of output file
        """

        t0 = time.time()

        #Note that a chunk is not a blob!!
        chunk_dim = self.__get_chunk_dimensions()
        blob_dim = self.__get_blob_dimensions(chunk_dim)
        n_blobs = self.container.calc_n_blobs(blob_dim)

        #Write header of .fil file
        n_bytes  = self.header['nbits'] / 8
        with open(filename_out, "w") as fileh:
            fileh.write(generate_sigproc_header(self))

        logger.info('Using %i n_blobs to write the data.'% n_blobs)
        for ii in range(0, n_blobs):
            logger.info('Reading %i of %i' % (ii + 1, n_blobs))

            bob = self.container.read_blob(blob_dim,n_blob=ii)

            #Write data of .fil file.
            with open(filename_out, "a") as fileh:
                j = bob
                if n_bytes == 4:
                    np.float32(j.ravel()).tofile(fileh)
                elif n_bytes == 2:
                    np.int16(j.ravel()).tofile(fileh)
                elif n_bytes == 1:
                    np.int8(j.ravel()).tofile(fileh)

        t1 = time.time()
        logger.info('Conversion time: %2.2fsec' % (t1- t0))

    def __write_to_fil_light(self, filename_out, *args, **kwargs):
        """ Write data to .fil file.

        Args:
            filename_out (str): Name of output file
        """

        t0 = time.time()

        n_bytes  = self.header['nbits'] / 8
        with open(filename_out, "w") as fileh:
            fileh.write(generate_sigproc_header(self))
            j = self.data
            if n_bytes == 4:
                np.float32(j.ravel()).tofile(fileh)
            elif n_bytes == 2:
                np.int16(j.ravel()).tofile(fileh)
            elif n_bytes == 1:
                np.int8(j.ravel()).tofile(fileh)

        t1 = time.time()
        logger.info('Conversion time: %2.2fsec' % (t1- t0))

    def write_to_hdf5(self, filename_out, *args, **kwargs):
        """ Write data to HDF5 file.
            It check the file size then decides how to write the file.

        Args:
            filename_out (str): Name of output file
        """

        if self.container.isheavy():
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

        #Note that a chunk is not a blob!!
        chunk_dim = self.__get_chunk_dimensions()
        blob_dim = self.__get_blob_dimensions(chunk_dim)
        n_blobs = self.container.calc_n_blobs(blob_dim)

        with h5py.File(filename_out, 'w') as h5:

            h5.attrs['CLASS'] = 'FILTERBANK'
            h5.attrs['VERSION'] = '1.0'

            if HAS_BITSHUFFLE:
                bs_compression = bitshuffle.h5.H5FILTER
                bs_compression_opts = (block_size, bitshuffle.h5.H5_COMPRESS_LZ4)
            else:
                bs_compression = None
                bs_compression_opts = None
                logger.warning("Warning: bitshuffle not found. No compression applied.")

            dset = h5.create_dataset('data',
                            shape=self.file_shape,
                            chunks=chunk_dim,
                            compression=bs_compression,
                            compression_opts=bs_compression_opts,
                            dtype=self.data.dtype)

            dset_mask = h5.create_dataset('mask',
                            shape=self.file_shape,
                            chunks=chunk_dim,
                            compression=bs_compression,
                            compression_opts=bs_compression_opts,
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

            if blob_dim[self.freq_axis] < self.n_channels_in_file:

                logger.info('Using %i n_blobs to write the data.'% n_blobs)
                for ii in range(0, n_blobs):
                    logger.info('Reading %i of %i' % (ii + 1, n_blobs))

                    bob = self.container.read_blob(blob_dim,n_blob=ii)

                    #-----
                    #Using channels instead of frequency.
                    c_start = self.container.c_start() + ii*blob_dim[self.freq_axis]
                    t_start = self.container.t_start + (c_start/self.n_channels_in_file)*blob_dim[self.time_axis]
                    t_stop = t_start + blob_dim[self.freq_axis]

                    # Reverse array if frequency axis is flipped
#                     if self.header['foff'] < 0:
#                         c_stop = self.n_channels_in_file - (c_start)%self.n_channels_in_file
#                         c_start = c_stop - blob_dim[self.freq_axis]
#                     else:
                    c_start = (c_start)%self.n_channels_in_file
                    c_stop = c_start + blob_dim[self.freq_axis]
                    #-----

                    logger.debug(t_start,t_stop,c_start,c_stop)

                    dset[t_start:t_stop,0,c_start:c_stop] = bob[:]

            else:

                logger.info('Using %i n_blobs to write the data.'% n_blobs)
                for ii in range(0, n_blobs):
                    logger.info('Reading %i of %i' % (ii + 1, n_blobs))

                    bob = self.container.read_blob(blob_dim,n_blob=ii)
                    t_start = self.container.t_start + ii*blob_dim[self.time_axis]
                    t_stop = min((ii+1)*blob_dim[self.time_axis],self.n_ints_in_file)

                    dset[t_start:t_stop] = bob[:]

        t1 = time.time()
        logger.info('Conversion time: %2.2fsec' % (t1- t0))

    def __write_to_hdf5_light(self, filename_out, *args, **kwargs):
        """ Write data to HDF5 file in one go.

        Args:
            filename_out (str): Name of output file
        """

        t0 = time.time()
        block_size = 0

        with h5py.File(filename_out, 'w') as h5:

            h5.attrs['CLASS'] = 'FILTERBANK'
            h5.attrs['VERSION'] = '1.0'

            if HAS_BITSHUFFLE:
                bs_compression = bitshuffle.h5.H5FILTER
                bs_compression_opts = (block_size, bitshuffle.h5.H5_COMPRESS_LZ4)
            else:
                bs_compression = None
                bs_compression_opts = None
                logger.warning("Warning: bitshuffle not found. No compression applied.")


            dset = h5.create_dataset('data',
                        data=self.data,
#                          compression='lzf')
                        compression=bs_compression,
                        compression_opts=bs_compression_opts)

            dset_mask = h5.create_dataset('mask',
                        shape=self.file_shape,
#                                 compression='lzf',
                        compression=bs_compression,
                        compression_opts=bs_compression_opts,
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

        t1 = time.time()
        logger.info('Conversion time: %2.2fsec' % (t1- t0))

    def __get_blob_dimensions(self, chunk_dim):
        """ Sets the blob dimmentions, trying to read around 256 MiB at a time.
            This is assuming a chunk is about 1 MiB.
        """

        if self.n_channels_in_file > chunk_dim[self.freq_axis]*MAX_BLOB_MB:
            freq_axis_size = self.n_channels_in_file
#             while freq_axis_size > chunk_dim[self.freq_axis]*MAX_BLOB_MB:
#                 freq_axis_size /= 2
            time_axis_size = 1
        else:
            freq_axis_size = self.n_channels_in_file
            time_axis_size = chunk_dim[self.time_axis] * MAX_BLOB_MB * chunk_dim[self.freq_axis] / freq_axis_size

        blob_dim = (time_axis_size, 1, freq_axis_size)

        return blob_dim

    def __get_chunk_dimensions(self):
        """ Sets the chunking dimmentions depending on the file type.
        """

        if '.0000.' in self.filename:
            logger.info('Detecting high frequency resolution data.')
            chunk_dim = (1,1,1048576) #1048576 is the number of channels in a coarse channel.
            return chunk_dim
        elif '.0001.' in self.filename:
            logger.info('Detecting high time resolution data.')
            chunk_dim = (2048,1,512) #512 is the total number of channels per single band (ie. blc00)
            return chunk_dim
        elif '.0002.' in self.filename:
            logger.info('Detecting intermediate frequency and time resolution data.')
            chunk_dim = (10,1,65536)  #65536 is the total number of channels per single band (ie. blc00)
#            chunk_dim = (1,1,65536/4)
            return chunk_dim
        else:
            logger.warning('File format not know. Will use minimum chunking. NOT OPTIMAL.')
            chunk_dim = (1,1,512)
            return chunk_dim

    def calc_n_coarse_chan(self):
        """ This makes an attempt to calculate the number of coarse channels in a given freq selection.
            It assumes for now that a single coarse channel is 2.9296875 MHz
        """

        n_coarse_chan = self.container.calc_n_coarse_chan()

        return n_coarse_chan

    def blank_dc(self, n_coarse_chan):
        """ Blank DC bins in coarse channels.

        Note: currently only works if entire waterfall file is read
        """

        n_chan = self.data.shape[-1]
        n_chan_per_coarse = n_chan / n_coarse_chan

        mid_chan = (n_chan_per_coarse / 2) - 1

        for ii in range(n_coarse_chan):
            ss = ii*n_chan_per_coarse
            self.data[..., ss+mid_chan] = np.median(self.data[..., ss+mid_chan+1:ss+mid_chan+10])


#EE Needs update
def cmd_tool(args=None):
    """ Command line tool for plotting and viewing info on waterfall files """

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Command line utility for reading and plotting waterfall files.")

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

    # Open blimpy data
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

    fil = Waterfall(filename, f_start=parse_args.f_start, f_stop=parse_args.f_stop,t_start=parse_args.t_start, t_stop=parse_args.t_stop,load_data=load_data)
    fil.info()

    #Check the size of selection.
    if fil.container.isheavy() or parse_args.to_hdf5 or parse_args.to_fil:
        info_only = True

    # And if we want to plot data, then plot data.

    if not info_only:
        print('')

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
            logger.info("Blanking DC bin")
            n_coarse_chan = fil.calc_n_coarse_chan()
            fil.blank_dc(n_coarse_chan)

        if parse_args.what_to_plot == "w":
            plt.figure("waterfall", figsize=(8, 6))
            fil.plot_waterfall(f_start=parse_args.f_start, f_stop=parse_args.f_stop)
        elif parse_args.what_to_plot == "s":
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
            fil.plot_time_series(f_start=parse_args.f_start, f_stop=parse_args.f_stop)
        elif parse_args.what_to_plot == "a":
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all')
        elif parse_args.what_to_plot == "ank":
            plt.figure("Multiple diagnostic plots", figsize=(12, 9),facecolor='white')
            fil.plot_all(logged=True, f_start=parse_args.f_start, f_stop=parse_args.f_stop, t='all',kutosis=False)

        if parse_args.plt_filename != '':
            plt.savefig(parse_args.plt_filename)

        if not parse_args.save_only:
            if 'DISPLAY' in os.environ.keys():
                plt.show()
            else:
                logger.warning("No $DISPLAY available.")


    else:

        if parse_args.to_hdf5 and parse_args.to_fil:
            raise Warning('Either provide to_hdf5 or to_fil, but not both.')

        if parse_args.to_hdf5:
            if not filename_out:
                filename_out = filename.replace('.fil','.h5')
            elif '.h5' not in filename_out:
                filename_out = filename_out.replace('.fil','')+'.h5'

            logger.info('Writing file : %s'%(filename_out))
            fil.write_to_hdf5(filename_out)
            logger.info('File written.')

        if parse_args.to_fil:
            if not filename_out:
                filename_out = filename.replace('.h5','.fil')
            elif '.fil' not in filename_out:
                filename_out = filename_out.replace('.h5','')+'.fil'

            logger.info('Writing file : %s'%(filename_out))
            fil.write_to_fil(filename_out)
            logger.info('File written.')


if __name__ == "__main__":
    cmd_tool()
