from .config import *
from . import plot_time_series, plot_kurtosis, plot_spectrum_min_max, plot_waterfall, plot_spectrum
from astropy import units as u


def plot_all(wf, t=0, f_start=None, f_stop=None, logged=False, if_id=0, kurtosis=True, **kwargs):
    """ Plot waterfall of data as well as spectrum; also, placeholder to make even more complicated plots in the future.

    Args:
        f_start (float): start frequency, in MHz
        f_stop (float): stop frequency, in MHz
        logged (bool): Plot in linear (False) or dB units (True),
        t (int): integration number to plot (0 -> len(data))
        logged (bool): Plot in linear (False) or dB units (True)
        if_id (int): IF identification (if multiple IF signals in file)
        kwargs: keyword args to be passed to matplotlib plot() and imshow()
    """
    if wf.header['nbits'] <= 2:
        logged = False

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.35, 0.5
    bottom, height = 0.45, 0.5
    width2, height2 = 0.1125, 0.15
    bottom2, left2 = bottom - height2 - .025, left - width2 - .02
    bottom3, left3 = bottom2 - height2 - .025, 0.075

    rect_waterfall = [left, bottom, width, height]
    rect_colorbar = [left + width, bottom, .025, height]
    rect_spectrum = [left, bottom2, width, height2]
    rect_min_max = [left, bottom3, width, height2]
    rect_timeseries = [left + width, bottom, width2, height]
    rect_kurtosis = [left3, bottom3, 0.25, height2]
    rect_header = [left3 - .05, bottom, 0.2, height]

    # --------
    #         axColorbar = plt.axes(rect_colorbar)
    #         print 'Ploting Colorbar'
    #         print plot_data.max()
    #         print plot_data.min()
    #
    #         plot_colorbar = range(plot_data.min(),plot_data.max(),int((plot_data.max()-plot_data.min())/plot_data.shape[0]))
    #         plot_colorbar = np.array([[plot_colorbar],[plot_colorbar]])
    #
    #         plt.imshow(plot_colorbar,aspect='auto', rasterized=True, interpolation='nearest',)

    #         axColorbar.xaxis.set_major_formatter(nullfmt)
    #         axColorbar.yaxis.set_major_formatter(nullfmt)

    #         heatmap = axColorbar.pcolor(plot_data, edgecolors = 'none', picker=True)
    #         plt.colorbar(heatmap, cax = axColorbar)
    # --------
    
    axMinMax = plt.axes(rect_min_max)
    print('Plotting Min Max')
    plot_spectrum_min_max(wf, logged=logged, f_start=f_start, f_stop=f_stop, t=t)
    plt.title('')
    axMinMax.yaxis.tick_right()
    axMinMax.yaxis.set_label_position("right")

    # --------
    axSpectrum = plt.axes(rect_spectrum,sharex=axMinMax)
    print('Plotting Spectrum')
    plot_spectrum(wf, logged=logged, f_start=f_start, f_stop=f_stop, t=t)
    plt.title('')
    axSpectrum.yaxis.tick_right()
    axSpectrum.yaxis.set_label_position("right")
    plt.xlabel('')
#        axSpectrum.xaxis.set_major_formatter(nullfmt)
    plt.setp(axSpectrum.get_xticklabels(), visible=False)

    # --------
    axWaterfall = plt.axes(rect_waterfall,sharex=axMinMax)
    print('Plotting Waterfall')
    plot_waterfall(wf, f_start=f_start, f_stop=f_stop, logged=logged, cb=False)
    plt.xlabel('')

    # no labels
#        axWaterfall.xaxis.set_major_formatter(nullfmt)
    plt.setp(axWaterfall.get_xticklabels(), visible=False)

    # --------
    axTimeseries = plt.axes(rect_timeseries)
    print('Plotting Timeseries')
    plot_time_series(wf, f_start=f_start, f_stop=f_stop, orientation='v')
    axTimeseries.yaxis.set_major_formatter(nullfmt)
#        axTimeseries.xaxis.set_major_formatter(nullfmt)

    # --------
    # Could exclude since it takes much longer to run than the other plots.
    if kurtosis:
        axKurtosis = plt.axes(rect_kurtosis)
        print('Plotting Kurtosis')
        plot_kurtosis(wf, f_start=f_start, f_stop=f_stop)


    # --------
    axHeader = plt.axes(rect_header)
    print('Plotting Header')
    # Generate nicer header
    telescopes = {0: 'Fake data',
                  1: 'Arecibo',
                  2: 'Ooty',
                  3: 'Nancay',
                  4: 'Parkes',
                  5: 'Jodrell',
                  6: 'GBT',
                  8: 'Effelsberg',
                  10: 'SRT',
                  64: 'MeerKAT',
                  65: 'KAT7'
                  }

    telescope = telescopes.get(wf.header['telescope_id'], wf.header['telescope_id'])

    plot_header = "%14s: %s\n" % ("TELESCOPE_ID", telescope)
    for key in ('SRC_RAJ', 'SRC_DEJ', 'TSTART', 'NCHANS', 'NBEAMS', 'NIFS', 'NBITS'):
        try:
            plot_header += "%14s: %s\n" % (key, wf.header[key.lower()])
        except KeyError:
            pass
    fch1 = "%6.6f MHz" % wf.header['fch1']

    foff = (wf.header['foff'] * 1e6 * u.Hz)
    if np.abs(foff) > 1e6 * u.Hz:
        foff = str(foff.to('MHz'))
    elif np.abs(foff) > 1e3 * u.Hz:
        foff = str(foff.to('kHz'))
    else:
        foff = str(foff.to('Hz'))

    plot_header += "%14s: %s\n" % ("FCH1", fch1)
    plot_header += "%14s: %s\n" % ("FOFF", foff)

    plt.text(0.05, .95, plot_header, ha='left', va='top', wrap=True)

    axHeader.set_facecolor('white')
    axHeader.xaxis.set_major_formatter(nullfmt)
    axHeader.yaxis.set_major_formatter(nullfmt)
