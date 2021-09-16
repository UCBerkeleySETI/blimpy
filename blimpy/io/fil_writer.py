"""
Procedures for writing to a Filterbank File
"""
import time
import numpy as np
from .sigproc import generate_sigproc_header


def write_to_fil(wf, filename_out):
    """ Write data to .fil file.
        It checks the file size then decides how to write the file.

    Args:
        wf : Waterfall object
        filename_out : str
            Name of output file
    """

    # For timing how long it takes to write a file.
    t0 = time.time()

    # Update header
    wf._update_header()

    if wf.container.isheavy():
        __write_to_fil_heavy(wf, filename_out)
    else:
        __write_to_fil_light(wf, filename_out)

    t1 = time.time()
    wf.logger.info('Conversion time: %2.2fsec' % (t1 - t0))


def __write_to_fil_heavy(wf, filename_out):
    """ Write data to .fil file.

    Args:
        wf : Waterfall object
        filename_out : str
            Name of output file
    """

    # Note that a chunk is not a blob!!
    chunk_dim = wf._get_chunk_dimensions()
    blob_dim = wf._get_blob_dimensions(chunk_dim)
    n_blobs = wf.container.calc_n_blobs(blob_dim)

    # Calculate number of bytes per data element
    n_bytes = wf.header['nbits'] / 8

    wf.logger.info("__write_to_fil_heavy: For {}, chunk_dim={}, blob_dim={}, n_blobs={}"
                   .format(filename_out, chunk_dim, blob_dim, n_blobs))

    with open(filename_out, "wb") as fileh:

        # Write header of .fil file
        fileh.write(generate_sigproc_header(wf))

        # For each blob
        for ii in range(0, n_blobs):

            wf.logger.info('__write_to_fil_heavy: Processing %i of %i' % (ii + 1, n_blobs))
            bob = wf.container.read_blob(blob_dim, n_blob=ii)

            # Write data of .fil file.
            if n_bytes == 4:
                np.float32(bob.ravel()).tofile(fileh)
            elif n_bytes == 2:
                np.int16(bob.ravel()).tofile(fileh)
            elif n_bytes == 1:
                np.int8(bob.ravel()).tofile(fileh)


def __write_to_fil_light(wf, filename_out):
    """ Write data to .fil file.

    Args:
        wf : Waterfall object
        filename_out : str
            Name of output file
    """

    wf.logger.info("__write_to_fil_light: Writing the spectra matrix for {} in one go."
                   .format(filename_out))
    n_bytes = wf.header['nbits'] / 8
    with open(filename_out, "wb") as fileh:
        fileh.write(generate_sigproc_header(wf))
        if n_bytes == 4:
            np.float32(wf.data.ravel()).tofile(fileh)
        elif n_bytes == 2:
            np.int16(wf.data.ravel()).tofile(fileh)
        elif n_bytes == 1:
            np.int8(wf.data.ravel()).tofile(fileh)
