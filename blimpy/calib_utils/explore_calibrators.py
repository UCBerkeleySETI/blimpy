from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = ascii.read('calibrators.txt')


def func_powerlaw(x,alph,c):
    """
    Power-law function
    """
    return c*x**alph

def get_freqs(MHz=False):
    """
    Return frequency array in calibrators.txt
    If MHz is true, then units of values are in MHz;
    Else the units are in GHz (default).
    """
    if MHz:
        return np.array(data['Freq']) * 1000
    return np.array(data['Freq'])

def cal_fluxes(source):
    """
    Return spectrum of a particular source (e.g. '3C295') in calibrators.txt
    """
    return np.array(data[source])

def errors():
    """
    Get errors in the data matrix.
    """
    return np.array(data['Err'])

def source_power_law_fit(source,minfreq,maxfreq):
    """
    Calculates optimized power-law parameters. Considers data for a particular source
    from calibrators.txt between minfreq and maxfreq (in MHz)
    """
    freqs = get_freqs(MHz=True) #In MHz
    fluxes = cal_fluxes(source)
    freqs_cut = freqs[np.where(np.logical_and(freqs>=minfreq, freqs<=maxfreq))]
    fluxes_cut = fluxes[np.where(np.logical_and(freqs>=minfreq, freqs<=maxfreq))]
    popt, dummy = curve_fit(func_powerlaw,freqs_cut,fluxes_cut)
    return popt

def plot_flux_comp(source,name=None,custom_minfreq=None,custom_maxfreq=None):
    """
    Plots the result of source_power_law_fit() along with the data from calibrators.txt
    for the source in question. Use 'name' to save the resulting plot.
    """
    freqs = get_freqs(MHz=True) #In MHz
    fluxes = cal_fluxes(source)
    errs = errors()
    if custom_minfreq is not None:
        minfreq = custom_minfreq
    else: minfreq = freqs.min()

    if custom_maxfreq is not None:
        maxfreq = custom_maxfreq
    else: maxfreq = freqs.max()

    alph,c = source_power_law_fit(source,minfreq,maxfreq)
    print(alph)

    freqvec = np.linspace(freqs.min(),freqs.max(),5000)
    model = func_powerlaw(freqvec,alph,c)

    plt.plot(freqvec,model,'r-',label='Power-law model')
    plt.errorbar(freqs,fluxes,yerr=fluxes*errs/100.0,marker='o',mfc='blue',markersize=2.5,capsize=2,capthick=0.5,fmt='go',label='Actual values')

    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux (Jy)')
    plt.yscale('log')
    plt.title('Predicted and Actual Spectrum of '+source)
    plt.legend()
    if name is not None:
        plt.savefig(name,dpi=250)
    plt.show()
