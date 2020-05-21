from .config import *
from ..utils import rebin, db

def plot_spectrum_min_max(wf, t=0, f_start=None, f_stop=None, logged=False, if_id=0, c=None, **kwargs):
    """ Plot frequency spectrum of a given file

    Args:
        logged (bool): Plot in linear (False) or dB units (True)
        if_id (int): IF identification (if multiple IF signals in file)
        c: color for line
        kwargs: keyword args to be passed to matplotlib plot()
    """
    ax = plt.gca()

    plot_f, plot_data = wf.grab_data(f_start, f_stop, if_id)

    # Using accending frequency for all plots.
    if wf.header['foff'] < 0:
        plot_data = plot_data[..., ::-1]  # Reverse data
        plot_f = plot_f[::-1]

    if logged:
        db_plot_data = db(plot_data[0])
        fig_max = np.nanmax(db_plot_data[db_plot_data != np.inf])
        fig_min = np.nanmin(db_plot_data[db_plot_data != -np.inf])
    else:
        fig_max = plot_data[0].max()
        fig_min = plot_data[0].min()

    print("averaging along time axis...")

    # Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
    if len(plot_data.shape) > 1:
        plot_max = plot_data.max(axis=0)
        plot_min = plot_data.min(axis=0)
        plot_data = plot_data.mean(axis=0)
    else:
        plot_max = plot_data.max()
        plot_min = plot_data.min()
        plot_data = plot_data.mean()

    # Rebin to max number of points
    dec_fac_x = 1
    MAX_PLT_POINTS = 8 * 64  # Low resoluition to see the difference.
    if plot_data.shape[0] > MAX_PLT_POINTS:
        dec_fac_x = int(plot_data.shape[0] / MAX_PLT_POINTS)

    plot_data = rebin(plot_data, dec_fac_x, 1)
    plot_min = rebin(plot_min, dec_fac_x, 1)
    plot_max = rebin(plot_max, dec_fac_x, 1)
    plot_f = rebin(plot_f, dec_fac_x, 1)

    if logged:
        plt.plot(plot_f, db(plot_data), "#333333", label='mean', **kwargs)
        plt.plot(plot_f, db(plot_max), "#e74c3c", label='max', **kwargs)
        plt.plot(plot_f, db(plot_min), '#3b5b92', label='min', **kwargs)
        plt.ylabel("Power [dB]")
    else:
        plt.plot(plot_f, plot_data, "#333333", label='mean', **kwargs)
        plt.plot(plot_f, plot_max, "#e74c3c", label='max', **kwargs)
        plt.plot(plot_f, plot_min, '#3b5b92', label='min', **kwargs)
        plt.ylabel("Power [counts]")
    plt.xlabel("Frequency [MHz]")
    plt.legend()

    try:
        plt.title(wf.header['source_name'])
    except KeyError:
        plt.title(wf.filename)

    plt.xlim(plot_f[0], plot_f[-1])
    if logged:
        try:
            plt.ylim(fig_min - 1, fig_max + 1)
        except ValueError:
            plt.ylim(-10, 20)
