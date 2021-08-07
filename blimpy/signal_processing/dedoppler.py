r""" De-doppler signal processing functions """


import numpy as np


def dedoppler_1(wf, drift_rate):
    """
    Simple de-doppler code for a Filterbank or HDF5 file.
    Parameters:
    ----------
    wf : object
        Blimpy Waterfall object, previously instantiated with a loaded data matrix.
    drift_rate : float
        Signal drift rate over time [Hz/s]
    """

    # Get the time sampling interval in seconda.
    tsamp = wf.header['tsamp']

    # Get the fine channel bandwidth in Hz.
    chan_bw = wf.header['foff'] * 1e6

    # Compute the number of numpy rolls to perform.
    n_roll = (drift_rate * tsamp) / chan_bw

    # For each time-row,
    #      roll all of the data power values in each fine channel frequency column
    #      given by -(n_roll * row number).
    for ii in range(wf.data.shape[0]):
        wf.data[ii][0][:] = np.roll(wf.data[ii][0][:], -int(n_roll * ii))
