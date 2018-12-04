import time
import h5py

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


def write_to_hdf5(wf, filename_out, *args, **kwargs):
    """ Write data to HDF5 file.
        It check the file size then decides how to write the file.

    Args:
        filename_out (str): Name of output file
    """

    #For timing how long it takes to write a file.
    t0 = time.time()

    #Update header
    wf._update_header()

    if wf.container.isheavy():
        __write_to_hdf5_heavy(wf, filename_out)
    else:
        __write_to_hdf5_light(wf, filename_out)

    t1 = time.time()
    wf.logger.info('Conversion time: %2.2fsec' % (t1- t0))


def __write_to_hdf5_heavy(wf, filename_out, *args, **kwargs):
    """ Write data to HDF5 file.

    Args:
        filename_out (str): Name of output file
    """

    block_size = 0

    #Note that a chunk is not a blob!!
    chunk_dim = wf._get_chunk_dimensions()
    blob_dim  = wf._get_blob_dimensions(chunk_dim)
    n_blobs   = wf.container.calc_n_blobs(blob_dim)

    with h5py.File(filename_out, 'w') as h5:

        h5.attrs[b'CLASS'] = b'FILTERBANK'
        h5.attrs[b'VERSION'] = b'1.0'

        if HAS_BITSHUFFLE:
            bs_compression = bitshuffle.h5.H5FILTER
            bs_compression_opts = (block_size, bitshuffle.h5.H5_COMPRESS_LZ4)
        else:
            bs_compression = None
            bs_compression_opts = None
            wf.logger.warning("Warning: bitshuffle not found. No compression applied.")

        dset = h5.create_dataset('data',
                                 shape=wf.selection_shape,
                                 chunks=chunk_dim,
                                 compression=bs_compression,
                                 compression_opts=bs_compression_opts,
                                 dtype=wf.data.dtype)

        dset_mask = h5.create_dataset('mask',
                                      shape=wf.selection_shape,
                                      chunks=chunk_dim,
                                      compression=bs_compression,
                                      compression_opts=bs_compression_opts,
                                      dtype='uint8')

        dset.dims[0].label = b"frequency"
        dset.dims[1].label = b"feed_id"
        dset.dims[2].label = b"time"

        dset_mask.dims[0].label = b"frequency"
        dset_mask.dims[1].label = b"feed_id"
        dset_mask.dims[2].label = b"time"

        # Copy over header information as attributes
        for key, value in wf.header.items():
            dset.attrs[key] = value

        if blob_dim[wf.freq_axis] < wf.selection_shape[wf.freq_axis]:

            wf.logger.info('Using %i n_blobs to write the data.'% n_blobs)
            for ii in range(0, n_blobs):
                wf.logger.info('Reading %i of %i' % (ii + 1, n_blobs))

                bob = wf.container.read_blob(blob_dim, n_blob=ii)

                #-----
                #Using channels instead of frequency.
                c_start = wf.container.chan_start_idx + ii * blob_dim[wf.freq_axis]
                t_start = wf.container.t_start + (c_start / wf.selection_shape[wf.freq_axis]) * blob_dim[wf.time_axis]
                t_stop = t_start + blob_dim[wf.time_axis]

                # Reverse array if frequency axis is flipped
#                     if self.header['foff'] < 0:
#                         c_stop = self.selection_shape[self.freq_axis] - (c_start)%self.selection_shape[self.freq_axis]
#                         c_start = c_stop - blob_dim[self.freq_axis]
#                     else:
                c_start = (c_start) % wf.selection_shape[wf.freq_axis]
                c_stop = c_start + blob_dim[wf.freq_axis]
                #-----

                wf.logger.debug(t_start,t_stop,c_start,c_stop)

                dset[t_start:t_stop,0,c_start:c_stop] = bob[:]

        else:

            wf.logger.info('Using %i n_blobs to write the data.'% n_blobs)
            for ii in range(0, n_blobs):
                wf.logger.info('Reading %i of %i' % (ii + 1, n_blobs))

                bob = wf.container.read_blob(blob_dim, n_blob=ii)
                t_start = wf.container.t_start + ii * blob_dim[wf.time_axis]

                #This prevents issues when the last blob is smaller than the others in time
                if (ii+1)*blob_dim[wf.time_axis] > wf.n_ints_in_file:
                    t_stop = wf.n_ints_in_file
                else:
                    t_stop = (ii+1)*blob_dim[wf.time_axis]

                dset[t_start:t_stop] = bob[:]


def __write_to_hdf5_light(wf, filename_out, *args, **kwargs):
    """ Write data to HDF5 file in one go.

    Args:
        filename_out (str): Name of output file
    """

    block_size = 0

    with h5py.File(filename_out, 'w') as h5:

        h5.attrs[b'CLASS']   = b'FILTERBANK'
        h5.attrs[b'VERSION'] = b'1.0'

        if HAS_BITSHUFFLE:
            bs_compression = bitshuffle.h5.H5FILTER
            bs_compression_opts = (block_size, bitshuffle.h5.H5_COMPRESS_LZ4)
        else:
            bs_compression = None
            bs_compression_opts = None
            wf.logger.warning("Warning: bitshuffle not found. No compression applied.")


        dset = h5.create_dataset('data',
                                 data=wf.data,
                                 #                          compression='lzf')
                                 compression=bs_compression,
                                 compression_opts=bs_compression_opts)

        dset_mask = h5.create_dataset('mask',
                                      shape=wf.file_shape,
                                      #                                 compression='lzf',
                                      compression=bs_compression,
                                      compression_opts=bs_compression_opts,
                                      dtype='uint8')

        dset.dims[0].label = b"frequency"
        dset.dims[1].label = b"feed_id"
        dset.dims[2].label = b"time"

        dset_mask.dims[0].label = b"frequency"
        dset_mask.dims[1].label = b"feed_id"
        dset_mask.dims[2].label = b"time"

        # Copy over header information as attributes
        for key, value in wf.header.items():
            dset.attrs[key] = value