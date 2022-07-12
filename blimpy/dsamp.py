"""
Downsample an input Filterbank file (.fil or .h5)
to an output .h5 Filterbank file.
"""


# External dependencies:
import sys
import time
from argparse import ArgumentParser
import numpy as np


# Logging set up:
import logging
LOGGER = logging.getLogger(__name__)
FMT = "%(name)-15s %(levelname)-8s %(message)s"
logging.basicConfig(format=FMT, stream=sys.stdout, level = logging.INFO)


# Blimpy functions required:
from blimpy import Waterfall
from blimpy.io.hdf_writer import __write_to_hdf5_heavy as write_to_h5


def downer(in_np_array, in_tsamp, group_size, out_dtype="float32"):
    """
    This is a downsample function.

    For every every group_size time samples of the input array,
         sum the element values into one total.
    The number of output samples = input array time dimension
         integer-divided by group_size.
    If the remainder of that division > 0,
         then the excess samples from the input array are dropped.

    Parameters
    ----------
    in_np_array : Numpy array
        Input data before summing.
    in_tsamp : float
        Input time sample size.
    group_size : int
        Group size for the purpose of summing 2 or more time samples.
    out_dtype : str
        Output data type.  Default is "float32".

    Returns
    -------
    Success:
        Downsampled data
        Output time sample size
        Output number of time integrations
    Failure: None, None, None.

    """
    # Get input array shape
    in_shape = in_np_array.shape
    if len(in_shape) != 3:
        LOGGER.error(f"Input array has {len(in_shape)} dimensions but 3 are required (time, nifs, fine-freqs) !!")
        return None, None, None
    if group_size < 2:
        LOGGER.error(f"Input group size ({group_size}) but it must be at least 2 !!")
        return None, None, None

    # Compute the number of sums.
    out_nints = np.floor_divide(in_shape[0], group_size)
    if out_nints < 1:
        LOGGER.error(f"Input group size ({group_size}) is larger than the time dimension of the input data ({in_shape[0]}) !!")
        return None, None, None
    LOGGER.info(f"Total input time samples to be dropped just before EOF = {in_shape[0] % group_size}")

    # Compute output time sample size.
    out_tsamp = in_tsamp * group_size

    # Initialise output array.
    out_np_array = np.zeros((out_nints, in_shape[1], in_shape[2]), dtype=out_dtype)

    # ii1 : time index that is bumped by group_size
    ii1 = 0

    # For each output row .....
    for mm in range(0, out_nints):

        # For each time row of the input array to sum for current output row .....
        for ii2 in range(ii1, ii1 + group_size):

            # For each polarisation in the row .....
            for jj in range(0, in_shape[1]):

                # For each find channel column in the polarisation .....
                for kk in range(0, in_shape[2]):

                    # Increment output element by an input element.
                    out_np_array[mm, jj, kk] += in_np_array[ii2, jj, kk]

        # Log progress.
        LOGGER.info(f"Completed {mm + 1} of {out_nints} output time samples.")

        # Point to the next group.
        ii1 += group_size

    # Done.  Return output array.
    return out_np_array, out_tsamp, out_nints


def make_h5_file(in_path, out_path, group_size):
    """
    1. Load input filterbank .fil or .h5 file.
    2. Call downer to perform down-sampling.
    3. Save result to the specified .h5 file.

    Args:
        in_path (str): Name of filterbank file to load
        out_path (str): Name of output filename. If not set, will default
                                    to same as input, but with .h5 instead of .fil
        group_size (int): Group size for the purpose of summing 2 or more time samples.
    """

    # Load input filterbank .fil or .h5 file.
    wf = Waterfall(in_path, load_data=True)

    # Down-sample input.
    t0 = time.time()
    out_data, out_tsamp, out_ntints = downer(wf.data, wf.header["tsamp"], group_size)
    if out_data is None:
        return 1
    LOGGER.info(f"Down-sampling time: {time.time() - t0 :f}s")
    LOGGER.info(f"Input data shape: {wf.data.shape}")

    # Write output file.
    wf.header["tsamp"] = out_tsamp
    wf.n_ints_in_file = out_ntints
    wf.selection_shape = (out_ntints, wf.header["nifs"], wf.n_channels_in_file)
    wf.file_shape =  wf.selection_shape
    wf.data = out_data
    LOGGER.info(f"Output data shape: {wf.data.shape}")
    t0 = time.time()
    write_to_h5(wf, out_path)
    LOGGER.info(f"Write-output time: {time.time() - t0 :f}s")
    return 0


def cmd_tool(args=None):
    """ Command line utility for downsampling a Filterbank file.
    """

    parser = ArgumentParser(description="Downsample an input Filterbank file (.fil or .h5) to an output .h5 Filterbank file.")
    parser.add_argument("in_path", type=str, help="Path of input Filterbank file (.fil or .h5)")
    parser.add_argument("out_path", type=str, help="Path of output Filterbank file (.h5 only)")
    parser.add_argument("-s", "--group_size", dest="group_size", type=int, required=True,
                        help="Group size for the purpose of summing 2 or more time samples.  Required.")

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    if args.group_size < 2:
        LOGGER.error(f"Input group size = {args.group_size} but it must be at least 2 !!")
        sys.exit(1)

    rc = make_h5_file(args.in_path,
                      args.out_path,
                      args.group_size)

    if rc != 0:
        sys.exit(rc)


if __name__ == "__main__":
    cmd_tool()
