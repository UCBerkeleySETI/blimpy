import os
import sys
import time
import h5py
import hdf5plugin
from blimpy import utils


def write_to_hdf5(wf, filename_out, f_scrunch=None, *args, **kwargs):
    """ Copy the header and the selected data matrix subset from the input file to the output HDF5 file.
        Check the heavy flag to decide how to write the file - light or heavy.

    Args:
        wf : Waterfall object
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
    """ Copy the header and the selected data matrix subset from the input file to the output HDF5 file, blob by blob.

    Args:
        wf : Waterfall object
        filename_out (str): Name of output file
        f_scrunch (int or None): Average (scrunch) N channels together
    """

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
    wf.logger.info("__write_to_hdf5_heavy: For {}, chunk_dim={}, blob_dim={}, n_blobs={}"
                   .format(filename_out, chunk_dim, blob_dim, n_blobs))

    # ===============================
    # Attempt to write the HDF5 file.
    # ===============================
   
    try:
        os.remove(filename_out)  # Try to pre-remove output .h5 file.
    except:
        pass
    
    try:
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
    
            # Copy over header information (Filterbank metadata) as HDF5 attributes
            for key, value in wf.header.items():
                dset.attrs[key] = value
    
            # Compare blob_dim along frequency axis
            # to the same of the data matrix selection.
            if blob_dim[wf.freq_axis] < wf.selection_shape[wf.freq_axis]:
                
                # Blob frequency axis is SHORTER THAN that of the data matrix selection.
    
                for ii in range(0, n_blobs):
                    wf.logger.info('__write_to_hdf5_heavy: Processing blob %i of %i' % (ii + 1, n_blobs))
    
                    # Read next data matrix blob from the input file.
                    bob = wf.container.read_blob(blob_dim, n_blob=ii)
    
                    # Using channel index values for c_start and c_stop.
                    c_start = wf.container.chan_start_idx + ii * blob_dim[wf.freq_axis]
                    t_start = wf.container.t_start + (c_start / wf.selection_shape[wf.freq_axis]) * blob_dim[wf.time_axis]
                    t_stop = t_start + blob_dim[wf.time_axis]  
                    c_start = (c_start) % wf.selection_shape[wf.freq_axis]
                    c_stop = c_start + blob_dim[wf.freq_axis]

                    # Frequency scrunching?    
                    if f_scrunch is not None:
                        c_start //= f_scrunch
                        c_stop  //= f_scrunch
                        bob = utils.rebin(bob, n_z=f_scrunch)
    
                    # Store blob in the output data.    
                    dset[t_start:t_stop, 0, c_start:c_stop] = bob[:]
                    wf.logger.debug('__write_to_hdf5_heavy t_start={}, t_stop={}, c_start={}, c_stop={}'
                                    .format(t_start, t_stop, c_start, c_stop))
    
            else:
                
                # Blob frequency axis is LONGER THAN or EQUAL TO that of the data matrix selection.
    
                for ii in range(0, n_blobs):
                    wf.logger.info('__write_to_hdf5_heavy: Processing blob %i of %i' % (ii + 1, n_blobs))
    
                    # Read next data matrix blob from the input file.
                    bob = wf.container.read_blob(blob_dim, n_blob=ii)
                    
                    # Compute start time.
                    t_start = wf.container.t_start + ii * blob_dim[wf.time_axis]
    
                    # Compute stop time.
                    if (ii + 1) * blob_dim[wf.time_axis] > wf.n_ints_in_file:
                        # Last blob: smaller than the others in the time dimension.
                        t_stop = wf.n_ints_in_file
                    else:
                        # All but the last blob: Full-sized blob in the time dimension.
                        t_stop = (ii + 1) * blob_dim[wf.time_axis]
    
                    # Frequency scrunching?    
                    if f_scrunch is not None:
                        bob = utils.rebin(bob, n_z=f_scrunch)

                    # Store blob in the output data.
                    # Note that the entire selected frequency dimension is used
                    # instead of the blob frequency dimension.
                    dset[t_start:t_stop] = bob[:]
                    wf.logger.debug('__write_to_hdf5_heavy t_start={}, t_stop={}'
                                    .format(t_start, t_stop))

        # =============================================
        # Success!  The .h5 file is written and closed.
        # =============================================

    except Exception as ex1: # Something went wrong!
        wf.logger.error("__write_to_hdf5_heavy: Writing the output HDF5 file {} failed!"
                       .format(filename_out))
        print(repr(ex1))
        try:
            if os.path.exists(filename_out):
                os.remove(filename_out) # scrap a potentially corrupted HDF5 file.
                # Removal succeeded.  Exit to the O/S with a nonzero exit code.
                wf.logger.info("__write_to_hdf5_heavy: Removal of partial HDF5 file {} succeeded."
                               .format(filename_out))
            sys.exit(86)
        except Exception as ex2:
            wf.logger.error("__write_to_hdf5_heavy: Removal of partial HDF5 file {} failed!"
                           .format(filename_out))
            print(repr(ex2))
            sys.exit(86)


def __write_to_hdf5_light(wf, filename_out, f_scrunch=None, *args, **kwargs):
    """ Copy the header and the selected data matrix subset from the input file to the output HDF5 file in one go.

    Args:
        wf : Waterfall object
        filename_out (str): Name of output file
        f_scrunch (int or None): Average (scrunch) N channels together
    """

    wf.logger.info("__write_to_hdf5_light: Writing the spectra matrix for {} without blobbing."
                   .format(filename_out))
    try:
        os.remove(filename_out)  # Try to pre-remove output .h5 file.
    except:
        pass # Ignore errors.
    
    with h5py.File(filename_out, 'w') as h5:

        h5.attrs['CLASS']   = 'FILTERBANK'
        h5.attrs['VERSION'] = '1.0'

        bs_compression = hdf5plugin.Bitshuffle(nelems=0, lz4=True)['compression']
        bs_compression_opts = hdf5plugin.Bitshuffle(nelems=0, lz4=True)['compression_opts']

        # Frequency scrunching?
        if f_scrunch is None:
            data_out = wf.data # No, just copy as-is from SIGPROC Filterbank (.fil) selection
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

        # Copy over header information (Filterbank metadata) as HDF5 attributes
        for key, value in wf.header.items():
            dset.attrs[key] = value
