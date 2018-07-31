import pylab as plt
from fluxcal import *
from stokescal import *
from blimpy import Waterfall
import numpy as np

def get_diff(dio_cross,**kwargs):
    '''
    Returns ON-OFF for all Stokes parameters given a cross_pols noise diode measurement
    '''
    #Get Stokes parameters, frequencies, and time sample length
    obs = Waterfall(dio_cross,max_load=150)
    freqs = obs.populate_freqs()
    tsamp = obs.header['tsamp']
    data = obs.data

    I,Q,U,V = get_stokes(data)

    #Fold noise diode data
    I_OFF,I_ON = foldcal(I,tsamp,**kwargs)
    Q_OFF,Q_ON = foldcal(Q,tsamp,**kwargs)
    U_OFF,U_ON = foldcal(U,tsamp,**kwargs)
    V_OFF,V_ON = foldcal(V,tsamp,**kwargs)

    #Do ON-OFF subtraction
    Idiff = I_ON-I_OFF
    Qdiff = Q_ON-Q_OFF
    Udiff = U_ON-U_OFF
    Vdiff = V_ON-V_OFF

    return Idiff,Qdiff,Udiff,Vdiff,freqs

def plot_Stokes_diode(dio_cross,diff=True,**kwargs):
    '''
    Plots the uncalibrate full stokes spectrum of the noise diode.
    Use diff=False to plot both ON and OFF, or diff=True for ON-OFF
    '''

    #If diff=True, get ON-OFF. If not get ON and OFF separately
    if diff==True:
        Idiff,Qdiff,Udiff,Vdiff,freqs = get_diff(dio_cross,**kwargs)
    else:
        obs = Waterfall(dio_cross,max_load=150)
        freqs = obs.populate_freqs()
        tsamp = obs.header['tsamp']
        data = obs.data
        I,Q,U,V = get_stokes(data)

        I_OFF,I_ON = foldcal(I,tsamp,**kwargs)
        Q_OFF,Q_ON = foldcal(Q,tsamp,**kwargs)
        U_OFF,U_ON = foldcal(U,tsamp,**kwargs)
        V_OFF,V_ON = foldcal(V,tsamp,**kwargs)

    #Plot spectra
    if diff==True:
        plt.plot(freqs,Idiff,'k-',label='I')
        plt.plot(freqs,Qdiff,'r-',label='Q')
        plt.plot(freqs,Udiff,'g-',label='U')
        plt.plot(freqs,Vdiff,'m-',label='V')
    else:
        plt.plot(freqs,I_ON,'k-',label='I ON')
        plt.plot(freqs,I_OFF,'k--',label='I OFF')
        plt.plot(freqs,Q_ON,'r-',label='Q ON')
        plt.plot(freqs,Q_OFF,'r--',label='Q OFF')
        plt.plot(freqs,U_ON,'g-',label='U ON')
        plt.plot(freqs,U_OFF,'g--',label='U OFF')
        plt.plot(freqs,V_ON,'m-',label='V ON')
        plt.plot(freqs,V_OFF,'m--',label='V OFF')

    plt.legend()
    plt.xlabel('Frequency (MHz)')
    plt.title('Uncalibrated Full Stokes Noise Diode Spectrum')
    plt.ylabel('Power (Counts)')
    plt.show()

def plot_calibrated_diode(dio_cross,chan_per_coarse=8,**kwargs):
    '''
    Plots the corrected noise diode spectrum for a given noise diode measurement
    after application of the inverse Mueller matrix for the electronics chain.
    '''
    #Get full stokes data for the ND observation
    obs = Waterfall(dio_cross,max_load=150)
    freqs = obs.populate_freqs()
    tsamp = obs.header['tsamp']
    data = obs.data
    I,Q,U,V = get_stokes(data)

    #Calculate Mueller Matrix variables for each coarse channel
    psis = phase_offsets(U,V,freqs,tsamp,chan_per_coarse,**kwargs)
    G = gain_offsets(I,Q,tsamp,chan_per_coarse,**kwargs)

    #Apply the Mueller matrix to original noise diode data and refold
    Icorr,Qcorr,Ucorr,Vcorr = apply_Mueller(I,Q,U,V,G,psis,chan_per_coarse)
    I_OFF,I_ON = foldcal(Icorr,tsamp,**kwargs)
    Q_OFF,Q_ON = foldcal(Qcorr,tsamp,**kwargs)
    U_OFF,U_ON = foldcal(Ucorr,tsamp,**kwargs)
    V_OFF,V_ON = foldcal(Vcorr,tsamp,**kwargs)

    #Plot new ON-OFF spectra
    plt.plot(freqs,I_ON-I_OFF,'k-',label='I')
    plt.plot(freqs,Q_ON-Q_OFF,'r-',label='Q')
    plt.plot(freqs,U_ON-U_OFF,'g-',label='U')
    plt.plot(freqs,V_ON-V_OFF,'m-',label='V')

    plt.legend()
    plt.xlabel('Frequency (MHz)')
    plt.title('Calibrated Full Stokes Noise Diode Spectrum')
    plt.ylabel('Power (Counts)')
    plt.show()

def plot_phase_offsets(dio_cross,chan_per_coarse=8,**kwargs):
    '''
    Plots the calculated phase offsets of each coarse channel along with
    the UV noise diode spectrum for comparison
    '''
    #Get ON-OFF ND spectra
    Idiff,Qdiff,Udiff,Vdiff,freqs = get_diff(dio_cross,**kwargs)
    obs = Waterfall(dio_cross,max_load=150)
    tsamp = obs.header['tsamp']
    data = obs.data
    I,Q,U,V = get_stokes(data)

    #Get phase offsets and convert to degrees
    coarse_psis = phase_offsets(U,V,freqs,tsamp,chan_per_coarse,**kwargs)
    coarse_freqs = convert_to_coarse(freqs,chan_per_coarse)
    coarse_degs = np.degrees(coarse_psis)

    #Plot phase offsets
    plt.subplot(211)
    plt.plot(coarse_freqs,coarse_degs,'ko',markersize=2,label='Coarse Channel $\psi$')
    plt.ylabel('Phase (Degrees)')
    plt.grid(True)
    plt.legend()
    plt.title('XY Phase Offsets')

    #Plot U and V spectra
    plt.subplot(212)
    plt.plot(freqs,Udiff,'g-',label='U')
    plt.plot(freqs,Vdiff,'m-',label='V')
    plt.legend()
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Power (Counts)')
    plt.grid(True)
    plt.show()


def plot_UVfit(dio_cross,chan_per_coarse=8,**kwargs):
    '''
    Plots phase offsets and sinusoidal fits of U and V for a noise diode measurement.
    NOTE: ONLY WORKS WHEN USING THE FIT METHOD OF FINDING PHASE OFFSETS.
    '''

    I,Q,U,V = get_stokes(dio_cross)
    tsamp = Waterfall(dio_cross,max_load=150).header['tsamp']

    Idiff,Qdiff,Udiff,Vdiff,freqs = get_diff(dio_cross,**kwargs)
    coarse_psis,Ufit,Vfit = phase_offsets(U,V,freqs,tsamp,chan_per_coarse,fit=True,**kwargs)
    coarse_freqs = convert_to_coarse(freqs,chan_per_coarse)

    coarse_degs = np.degrees(coarse_psis)
    coarse_degs = coarse_degs % 360

    plt.subplot(211)
    plt.plot(coarse_freqs,coarse_degs,'ko',label='Fit Phase')
    plt.ylabel('Phase (Degrees)')
    plt.grid(True)
    plt.legend()
    plt.title('UV Noise Diode Fits')
    plt.yticks(np.arange(0,360,step=45))

    plt.subplot(212)
    plt.plot(freqs,Udiff,'g--',label='U')
    plt.plot(freqs,Vdiff,'m--',label='V')
    plt.plot(freqs,Ufit,'k-',label='U fit')
    plt.plot(freqs,Vfit,'k-',label='V fit')
    plt.legend()
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Power (Counts')

    plt.show()

def plot_diode_fold(dio_cross,**kwargs):
    '''
    Plots the calculated average power and time sampling of ON (red) and
    OFF (blue) for a noise diode measurement over the observation time series
    '''
    #Get full stokes data of ND measurement
    obs = Waterfall(dio_cross,max_load=150)
    tsamp = obs.header['tsamp']
    data = obs.data
    I,Q,U,V = get_stokes(data)

    #Calculate time series, OFF and ON averages, and time samples for each
    tseries = np.squeeze(np.mean(I,axis=2))
    I_OFF,I_ON,OFFints,ONints = foldcal(I,tsamp,inds=True,**kwargs)
    stop = ONints[-1,1]

    #Plot time series and calculated average for ON and OFF
    plt.plot(tseries[0:stop],'k-',label='Total Power')
    for i in ONints:
        plt.plot(np.arange(i[0],i[1]),np.full((i[1]-i[0]),np.mean(I_ON)),'r-')
    for i in OFFints:
        plt.plot(np.arange(i[0],i[1]),np.full((i[1]-i[0]),np.mean(I_OFF)),'b-')

    plt.legend()
    plt.xlabel('Time Sample Number')
    plt.ylabel('Power (Counts)')
    plt.title('Noise Diode Fold')
    plt.show()


