from .config import *
from ..utils import rebin, db
from .plot_utils import calc_extent


def plot_waterfall(wf, f_start=None, f_stop=None, if_id=0, logged=True, cb=True, MJD_time=False, **kwargs):
    """ Plot waterfall of data

    Args:
        f_start (float): start frequency, in MHz
        f_stop (float): stop frequency, in MHz
        logged (bool): Plot in linear (False) or dB units (True),
        cb (bool): for plotting the colorbar
        kwargs: keyword args to be passed to matplotlib imshow()
    """
    plot_f, plot_data = wf.grab_data(f_start, f_stop, if_id)
    
    # imshow does not support int8, so convert to floating point
    plot_data = plot_data.astype('float32')
    
    # Using accending frequency for all plots.
    if wf.header['foff'] < 0:
        plot_data = plot_data[..., ::-1]  # Reverse data
        plot_f = plot_f[::-1]

    if logged:
        plot_data = db(plot_data)

    # Make sure waterfall plot is under 4k*4k
    dec_fac_x, dec_fac_y = 1, 1
    if plot_data.shape[0] > MAX_IMSHOW_POINTS[0]:
        dec_fac_x = int(plot_data.shape[0] / MAX_IMSHOW_POINTS[0])

    if plot_data.shape[1] > MAX_IMSHOW_POINTS[1]:
        dec_fac_y = int(plot_data.shape[1] / MAX_IMSHOW_POINTS[1])

    plot_data = rebin(plot_data, dec_fac_x, dec_fac_y)

    try:
        plt.title(wf.header['source_name'])
    except KeyError:
        plt.title(wf.filename)

    extent = calc_extent(wf, plot_f=plot_f, plot_t=wf.timestamps, MJD_time=MJD_time)

    plt.imshow(plot_data,
               aspect='auto',
               origin='lower',
               rasterized=True,
               interpolation='nearest',
               extent=extent,
               cmap='viridis',
               **kwargs
               )
    if cb:
        plt.colorbar()
    plt.xlabel("Frequency [MHz]")
    if MJD_time:
        plt.ylabel("Time [MJD]")
    else:
        plt.ylabel("Time [s]")
