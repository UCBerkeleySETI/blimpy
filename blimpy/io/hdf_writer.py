import time
import h5py
import hdf5plugin
from blimpy import utils


def write_to_hdf5(wf, filename_out, f_scrunch=None, *args, **kwargs):
    """ Write data to HDF5 file.
        It check the file size then decides how to write the file.

    Args:
        filename_out (str): Name of output file
        f_scrunch (int or None): Average (scrunch) N channels together
    """

    #For timing how long it takes to write a file.
    t0 = time.time()

    #Update header
    wf._update_header()

    if wf.container.isheavy():
        __write_to_hdf5_heavy(wf, filename_out, f_scrunch=f_scrunch, *args, **kwargs)
    else:
        __write_to_hdf5_light(wf, filename_out, f_scrunch=f_scrunch, *args, **kwargs)

    t1 = time.time()
    wf.logger.info('Conversion time: %2.2fsec' % (t1- t0))


def __write_to_hdf5_heavy(wf, filename_out, f_scrunch=None, *args, **kwargs):
    """ Write data to HDF5 file.

    Args:
        filename_out (str): Name of output file
        f_scrunch (int or None): Average (scrunch) N channels together
    """

    block_size = 0

    # Note that a chunk is not a blob!!
    # chunk_dim = wf._get_chunk_dimensions() <-- seems intended for raw to fil
    # And, chunk dimensions should not exceed the Waterfall selection shape dimensions.
    chunk_list = list(wf._get_chunk_dimensions())
    for ix in range(0, len(chunk_list)):
        if chunk_list[ix] > wf.selection_shape[ix]:
            chunk_list[ix] = wf.selection_shape[ix]
    chunk_dim = tuple(chunk_list)
    blob_dim  = wf._get_blob_dimensions(chunk_dim)
    n_blobs   = wf.container.calc_n_blobs(blob_dim)

    with h5py.File(filename_out, 'w') as h5:

        h5.attrs['CLASS'] = 'FILTERBANK'
        h5.attrs['VERSION'] = '1.0'

        bs_compression = hdf5plugin.Bitshuffle(nelems=0, lz4=True)['compression']
        bs_compression_opts = hdf5plugin.Bitshuffle(nelems=0, lz4=True)['compression_opts']

        dout_shape     = list(wf.selection_shape)    # Make sure not a tuple
        dout_chunk_dim = list(chunk_dim)

        if f_scrunch is not None:
            dout_shape[-1] //= f_scrunch
            dout_chunk_dim[-1] //= f_scrunch
            wf.header['foff'] *= f_scrunch

        dset = h5.create_dataset('data',
                                 shape=tuple(dout_shape),
                                 chunks=tuple(dout_chunk_dim),
                                 compression=bs_compression,
                                 compression_opts=bs_compression_opts,
                                 dtype=wf.data.dtype)

        dset_mask = h5.create_dataset('mask',
                                      shape=tuple(dout_shape),
                                      chunks=tuple(dout_chunk_dim),
                                      compression=bs_compression,
                                      compression_opts=bs_compression_opts,
                                      dtype='uint8')

        dset.dims[2].label = b"frequency"
        dset.dims[1].label = b"feed_id"
        dset.dims[0].label = b"time"

        dset_mask.dims[2].label = b"frequency"
        dset_mask.dims[1].label = b"feed_id"
        dset_mask.dims[0].label = b"time"

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

                if f_scrunch is not None:
                    c_start //= f_scrunch
                    c_stop  //= f_scrunch
                    bob = utils.rebin(bob, n_z=f_scrunch)

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

                if f_scrunch is not None:
                    bob = utils.rebin(bob, n_z=f_scrunch)

                dset[t_start:t_stop] = bob[:]


def __write_to_hdf5_light(wf, filename_out, f_scrunch=None, *args, **kwargs):
    """ Write data to HDF5 file in one go.

    Args:
        filename_out (str): Name of output file
        f_scrunch (int or None): Average (scrunch) N channels together
    """

    block_size = 0

    with h5py.File(filename_out, 'w') as h5:

        h5.attrs['CLASS']   = 'FILTERBANK'
        h5.attrs['VERSION'] = '1.0'

        bs_compression = hdf5plugin.Bitshuffle(nelems=0, lz4=True)['compression']
        bs_compression_opts = hdf5plugin.Bitshuffle(nelems=0, lz4=True)['compression_opts']

        if f_scrunch is None:
            data_out = wf.data
        else:
            wf.logger.info('Frequency scrunching by %i' % f_scrunch)
            data_out = utils.rebin(wf.data, n_z=f_scrunch)
            wf.header['foff'] *= f_scrunch

        dset = h5.create_dataset('data',
                                 data=data_out,
                                 compression=bs_compression,
                                 compression_opts=bs_compression_opts)

        dset_mask = h5.create_dataset('mask',
                                      shape=data_out.shape,
                                      compression=bs_compression,
                                      compression_opts=bs_compression_opts,
                                      dtype='uint8')

        dset.dims[2].label = b"frequency"
        dset.dims[1].label = b"feed_id"
        dset.dims[0].label = b"time"

        dset_mask.dims[2].label = b"frequency"
        dset_mask.dims[1].label = b"feed_id"
        dset_mask.dims[0].label = b"time"

        # Copy over header information as attributes
        for key, value in wf.header.items():
            dset.attrs[key] = value

