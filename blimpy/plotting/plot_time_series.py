from .config import *
from ..utils import rebin, db
from .plot_utils import calc_extent

def plot_time_series(wf, f_start=None, f_stop=None, if_id=0, logged=True, orientation='h', MJD_time=False, **kwargs):
    """ Plot the time series.

     Args:
        f_start (float): start frequency, in MHz
        f_stop (float): stop frequency, in MHz
        logged (bool): Plot in linear (False) or dB units (True),
        kwargs: keyword args to be passed to matplotlib imshow()
    """

    ax = plt.gca()
    plot_f, plot_data = wf.grab_data(f_start, f_stop, if_id)

    # Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
    if len(plot_data.shape) > 1:
        plot_data = np.nanmean(plot_data, axis=1)
    else:
        plot_data = np.nanmean(plot_data)

    if logged and wf.header['nbits'] >= 8:
        plot_data = db(plot_data)

    # Make proper time axis for plotting (but only for plotting!). Note that this makes the values inclusive.
    extent = calc_extent(wf, plot_f=plot_f, plot_t=wf.timestamps, MJD_time=MJD_time)
    plot_t = np.linspace(extent[2], extent[3], len(wf.timestamps))

    if MJD_time:
        tlabel = "Time [MJD]"
    else:
        tlabel = "Time [s]"

    if logged:
        plabel = "Power [dB]"
    else:
        plabel = "Power [counts]"

    # Reverse oder if vertical orientation.
    if 'v' in orientation:
        plt.plot(plot_data, plot_t, **kwargs)
        plt.xlabel(plabel)

    else:
        plt.plot(plot_t, plot_data, **kwargs)
        plt.xlabel(tlabel)
        plt.ylabel(plabel)

    ax.autoscale(axis='both', tight=True)
