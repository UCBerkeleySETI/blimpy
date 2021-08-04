import numpy as np


class WaterfallBlob(np.ndarray):
    """
    Subclassing the waterfall numpy array to add functionalities
    like dedoppler. 

    Parameters
    ----------
    waterfall : np.ndarray
        waterfall numpy array

    header : dict
        header dictionary

    """
    def __new__(cls, waterfall, header):
        obj = waterfall.view(cls)
        obj.header = header
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        if hasattr(obj, "header"):
            self.header = obj.header

    def dedoppler(self, drift_rate):
        """ 
        Simple de-doppler code for a filterbank 

        Parameters:
        ----------
        drift_rate : float
            drift_rate [Hz/s]

        """
        tsamp   = self.header['tsamp']
        chan_bw = self.header['foff'] * 1e6 #foff is in MHz
        n_roll = (drift_rate * tsamp) / chan_bw
        for ii in range(self.shape[0]):
            self[ii] = np.roll(self[ii], -int(n_roll*ii))

