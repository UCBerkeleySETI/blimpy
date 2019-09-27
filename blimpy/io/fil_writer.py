from .sigproc import *
import time
import numpy as np


def write_to_fil(wf, filename_out, *args, **kwargs):
    """ Write data to .fil file.
        It check the file size then decides how to write the file.

    Args:
        filename_out (str): Name of output file
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


def __write_to_fil_heavy(wf, filename_out, *args, **kwargs):
    """ Write data to .fil file.

    Args:
        filename_out (str): Name of output file
    """

    # Note that a chunk is not a blob!!
    chunk_dim = wf._get_chunk_dimensions()
    blob_dim = wf._get_blob_dimensions(chunk_dim)
    n_blobs = wf.container.calc_n_blobs(blob_dim)

    # Write header of .fil file
    n_bytes = wf.header[b'nbits'] / 8
    with open(filename_out, "wb") as fileh:
        fileh.write(generate_sigproc_header(wf))  # generate_sigproc_header comes from sigproc.py

    wf.logger.info('Using %i n_blobs to write the data.' % n_blobs)
    for ii in range(0, n_blobs):
        wf.logger.info('Reading %i of %i' % (ii + 1, n_blobs))

        bob = wf.container.read_blob(blob_dim, n_blob=ii)

        # Write data of .fil file.
        with open(filename_out, "a") as fileh:
            j = bob
            if n_bytes == 4:
                np.float32(j.ravel()).tofile(fileh)
            elif n_bytes == 2:
                np.int16(j.ravel()).tofile(fileh)
            elif n_bytes == 1:
                np.int8(j.ravel()).tofile(fileh)


def __write_to_fil_light(wf, filename_out, *args, **kwargs):
    """ Write data to .fil file.

    Args:
        filename_out (str): Name of output file
    """

    n_bytes = wf.header[b'nbits'] / 8
    with open(filename_out, "wb") as fileh:
        fileh.write(generate_sigproc_header(wf))  # generate_sigproc_header comes from sigproc.py
        j = wf.data
        if n_bytes == 4:
            np.float32(j.ravel()).tofile(fileh)
        elif n_bytes == 2:
            np.int16(j.ravel()).tofile(fileh)
        elif n_bytes == 1:
            np.int8(j.ravel()).tofile(fileh)