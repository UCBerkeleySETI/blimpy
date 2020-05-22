from .config import *
import scipy.stats


def plot_kurtosis(wf, f_start=None, f_stop=None, if_id=0, **kwargs):
    """ Plot kurtosis

     Args:
        f_start (float): start frequency, in MHz
        f_stop (float): stop frequency, in MHz
        kwargs: keyword args to be passed to matplotlib imshow()
    """
    ax = plt.gca()

    plot_f, plot_data = wf.grab_data(f_start, f_stop, if_id)

    # Using accending frequency for all plots.
    if wf.header['foff'] < 0:
        plot_data = plot_data[..., ::-1]  # Reverse data
        plot_f = plot_f[::-1]

    try:
        pltdata = scipy.stats.kurtosis(plot_data, axis=0, nan_policy='omit')
    except:
        pltdata = plot_data * 0.0

    plt.plot(plot_f, pltdata, **kwargs)
    plt.ylabel("Kurtosis")
    plt.xlabel("Frequency [MHz]")

    plt.xlim(plot_f[0], plot_f[-1])
