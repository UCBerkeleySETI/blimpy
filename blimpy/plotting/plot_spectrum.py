from .config import *
from ..utils import rebin, db


def plot_spectrum(wf, t=0, f_start=None, f_stop=None, logged=False, if_id=0, c=None, **kwargs):
    """ Plot frequency spectrum of a given file

    Args:
        t (int): integration number to plot (0 -> len(data))
        logged (bool): Plot in linear (False) or dB units (True)
        if_id (int): IF identification (if multiple IF signals in file)
        c: color for line
        kwargs: keyword args to be passed to matplotlib plot()
    """
    if wf.header['nbits'] <= 2:
        logged = False
        t = 'all'
    ax = plt.gca()

    plot_f, plot_data = wf.grab_data(f_start, f_stop, if_id)

    # Using accending frequency for all plots.
    if wf.header['foff'] < 0:
        plot_data = plot_data[..., ::-1]  # Reverse data
        plot_f = plot_f[::-1]

    if isinstance(t, int):
        print("extracting integration %i..." % t)
        plot_data = plot_data[t]
    elif t == 'all':
        print("averaging along time axis...")
        # Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
        if len(plot_data.shape) > 1:
            plot_data = plot_data.mean(axis=0)
        else:
            plot_data = plot_data.mean()
    else:
        raise RuntimeError("Unknown integration %s" % t)

    # Rebin to max number of points
    dec_fac_x = 1
    if plot_data.shape[0] > MAX_PLT_POINTS:
        dec_fac_x = int(plot_data.shape[0] / MAX_PLT_POINTS)

    plot_data = rebin(plot_data, dec_fac_x, 1)
    plot_f = rebin(plot_f, dec_fac_x, 1)

    if not c:
        kwargs['c'] = '#333333'

    if logged:
        plt.plot(plot_f, db(plot_data), label='Stokes I', **kwargs)
        plt.ylabel("Power [dB]")
    else:

        plt.plot(plot_f, plot_data, label='Stokes I', **kwargs)
        plt.ylabel("Power [counts]")
    plt.xlabel("Frequency [MHz]")
    plt.legend()

    try:
        plt.title(wf.header['source_name'])
    except KeyError:
        plt.title(wf.filename)

    plt.xlim(plot_f[0], plot_f[-1])
