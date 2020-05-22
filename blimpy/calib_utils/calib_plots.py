import pylab as plt
from .fluxcal import *
from .stokescal import *
from blimpy import Waterfall
import numpy as np

def get_diff(dio_cross,feedtype,**kwargs):
    """
    Returns ON-OFF for all Stokes parameters given a cross_pols noise diode measurement
    """
    #Get Stokes parameters, frequencies, and time sample length
    obs = Waterfall(dio_cross,max_load=150)
    freqs = obs.populate_freqs()
    tsamp = obs.header['tsamp']
    data = obs.data
    obs = None

    I,Q,U,V = get_stokes(data,feedtype)

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

def plot_Stokes_diode(dio_cross,diff=True,feedtype='l',**kwargs):
    """
    Plots the uncalibrated full stokes spectrum of the noise diode.
    Use diff=False to plot both ON and OFF, or diff=True for ON-OFF
    """

    #If diff=True, get ON-OFF. If not get ON and OFF separately
    if diff==True:
        Idiff,Qdiff,Udiff,Vdiff,freqs = get_diff(dio_cross,feedtype,**kwargs)
    else:
        obs = Waterfall(dio_cross,max_load=150)
        freqs = obs.populate_freqs()
        tsamp = obs.header['tsamp']
        data = obs.data
        I,Q,U,V = get_stokes(data,feedtype)

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


def plot_calibrated_diode(dio_cross,chan_per_coarse=8,feedtype='l',**kwargs):
    """
    Plots the corrected noise diode spectrum for a given noise diode measurement
    after application of the inverse Mueller matrix for the electronics chain.
    """
    #Get full stokes data for the ND observation
    obs = Waterfall(dio_cross,max_load=150)
    freqs = obs.populate_freqs()
    tsamp = obs.header['tsamp']
    data = obs.data
    obs = None
    I,Q,U,V = get_stokes(data,feedtype)
    data = None

    #Calculate Mueller Matrix variables for each coarse channel
    psis = phase_offsets(I,Q,U,V,tsamp,chan_per_coarse,feedtype,**kwargs)
    G = gain_offsets(I,Q,U,V,tsamp,chan_per_coarse,feedtype,**kwargs)

    #Apply the Mueller matrix to original noise diode data and refold
    I,Q,U,V = apply_Mueller(I,Q,U,V,G,psis,chan_per_coarse,feedtype)
    I_OFF,I_ON = foldcal(I,tsamp,**kwargs)
    Q_OFF,Q_ON = foldcal(Q,tsamp,**kwargs)
    U_OFF,U_ON = foldcal(U,tsamp,**kwargs)
    V_OFF,V_ON = foldcal(V,tsamp,**kwargs)

    #Delete data arrays for space
    I = None
    Q = None
    U = None
    V = None

    #Plot new ON-OFF spectra
    plt.plot(freqs,I_ON-I_OFF,'k-',label='I')
    plt.plot(freqs,Q_ON-Q_OFF,'r-',label='Q')
    plt.plot(freqs,U_ON-U_OFF,'g-',label='U')
    plt.plot(freqs,V_ON-V_OFF,'m-',label='V')

    plt.legend()
    plt.xlabel('Frequency (MHz)')
    plt.title('Calibrated Full Stokes Noise Diode Spectrum')
    plt.ylabel('Power (Counts)')


def plot_phase_offsets(dio_cross,chan_per_coarse=8,feedtype='l',ax1=None,ax2=None,legend=True,**kwargs):
    """
    Plots the calculated phase offsets of each coarse channel along with
    the UV (or QU) noise diode spectrum for comparison
    """
    #Get ON-OFF ND spectra
    Idiff,Qdiff,Udiff,Vdiff,freqs = get_diff(dio_cross,feedtype,**kwargs)
    obs = Waterfall(dio_cross,max_load=150)
    tsamp = obs.header['tsamp']
    data = obs.data
    obs = None
    I,Q,U,V = get_stokes(data,feedtype)

    #Get phase offsets and convert to degrees
    coarse_psis = phase_offsets(I,Q,U,V,tsamp,chan_per_coarse,feedtype,**kwargs)
    coarse_freqs = convert_to_coarse(freqs,chan_per_coarse)
    coarse_degs = np.degrees(coarse_psis)

    #Plot phase offsets
    if ax2==None:
        plt.subplot(211)
    else:
        axPsi = plt.axes(ax2)
        plt.setp(axPsi.get_xticklabels(),visible=False)
    plt.plot(coarse_freqs,coarse_degs,'ko',markersize=2,label='Coarse Channel $\psi$')
    plt.ylabel('Degrees')
    plt.grid(True)
    plt.title('Phase Offsets')
    if legend==True:
        plt.legend()

    #Plot U and V spectra
    if ax1==None:
        plt.subplot(212)
    else:
        axUV = plt.axes(ax1)
    plt.plot(freqs,Udiff,'g-',label='U')
    if feedtype=='l':
        plt.plot(freqs,Vdiff,'m-',label='V')
    if feedtype=='c':
        plt.plot(freqs,Qdiff,'r-',label='Q')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Power (Counts)')
    plt.grid(True)
    if legend==True:
        plt.legend()

def plot_gain_offsets(dio_cross,dio_chan_per_coarse=8,feedtype='l',ax1=None,ax2=None,legend=True,**kwargs):
    """
    Plots the calculated gain offsets of each coarse channel along with
    the time averaged power spectra of the X and Y feeds
    """
    #Get ON-OFF ND spectra
    Idiff,Qdiff,Udiff,Vdiff,freqs = get_diff(dio_cross,feedtype,**kwargs)
    obs = Waterfall(dio_cross,max_load=150)
    tsamp = obs.header['tsamp']
    data = obs.data
    obs = None
    I,Q,U,V = get_stokes(data,feedtype)

    #Get phase offsets and convert to degrees
    coarse_G = gain_offsets(I,Q,U,V,tsamp,dio_chan_per_coarse,feedtype,**kwargs)
    coarse_freqs = convert_to_coarse(freqs,dio_chan_per_coarse)

    #Get X and Y spectra for the noise diode ON and OFF
    #If using circular feeds these correspond to LL and RR
    XX_OFF,XX_ON = foldcal(np.expand_dims(data[:,0,:],axis=1),tsamp,**kwargs)
    YY_OFF,YY_ON = foldcal(np.expand_dims(data[:,1,:],axis=1),tsamp,**kwargs)

    if ax1==None:
        plt.subplot(211)
    else:
        axG = plt.axes(ax1)
        plt.setp(axG.get_xticklabels(),visible=False)
    plt.plot(coarse_freqs,coarse_G,'ko',markersize=2)
    plt.ylabel(r'$\frac{\Delta G}{2}$',rotation=90)
    if feedtype=='l':
        plt.title('XY Gain Difference')
    if feedtype=='c':
        plt.title('LR Gain Difference')
    plt.grid(True)

    if ax2==None:
        plt.subplot(212)
    else:
        axXY = plt.axes(ax2,sharex=axG)
    if feedtype=='l':
        plt.plot(freqs,XX_OFF,'b-',label='XX')
        plt.plot(freqs,YY_OFF,'r-',label='YY')
    if feedtype=='c':
        plt.plot(freqs,XX_OFF,'b-',label='LL')
        plt.plot(freqs,YY_OFF,'r-',label='RR')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Power (Counts)')
    if legend==True:
        plt.legend()

def plot_diode_fold(dio_cross,bothfeeds=True,feedtype='l',min_samp=-500,max_samp=7000,legend=True,**kwargs):
    """
    Plots the calculated average power and time sampling of ON (red) and
    OFF (blue) for a noise diode measurement over the observation time series
    """
    #Get full stokes data of ND measurement
    obs = Waterfall(dio_cross,max_load=150)
    tsamp = obs.header['tsamp']
    data = obs.data
    obs = None
    I,Q,U,V = get_stokes(data,feedtype)

    #Calculate time series, OFF and ON averages, and time samples for each
    tseriesI = np.squeeze(np.mean(I,axis=2))
    I_OFF,I_ON,OFFints,ONints = foldcal(I,tsamp,inds=True,**kwargs)
    if bothfeeds==True:
        if feedtype=='l':
            tseriesQ = np.squeeze(np.mean(Q,axis=2))
            tseriesX = (tseriesI+tseriesQ)/2
            tseriesY = (tseriesI-tseriesQ)/2
        if feedtype=='c':
            tseriesV = np.squeeze(np.mean(V,axis=2))
            tseriesR = (tseriesI+tseriesV)/2
            tseriesL = (tseriesI-tseriesV)/2

    stop = ONints[-1,1]

    #Plot time series and calculated average for ON and OFF
    if bothfeeds==False:
        plt.plot(tseriesI[0:stop],'k-',label='Total Power')
        for i in ONints:
            plt.plot(np.arange(i[0],i[1]),np.full((i[1]-i[0]),np.mean(I_ON)),'r-')
        for i in OFFints:
            plt.plot(np.arange(i[0],i[1]),np.full((i[1]-i[0]),np.mean(I_OFF)),'b-')
    else:
        if feedtype=='l':
            diff = np.mean(tseriesX)-np.mean(tseriesY)
            plt.plot(tseriesX[0:stop],'b-',label='XX')
            plt.plot(tseriesY[0:stop]+diff,'r-',label='YY (shifted)')
        if feedtype=='c':
            diff = np.mean(tseriesL)-np.mean(tseriesR)
            plt.plot(tseriesL[0:stop],'b-',label='LL')
            plt.plot(tseriesR[0:stop]+diff,'r-',label='RR (shifted)')

    #Calculate plotting limits
    if bothfeeds==False:
        lowlim = np.mean(I_OFF)-(np.mean(I_ON)-np.mean(I_OFF))/2
        hilim = np.mean(I_ON)+(np.mean(I_ON)-np.mean(I_OFF))/2
        plt.ylim((lowlim,hilim))

    plt.xlim((min_samp,max_samp))
    plt.xlabel('Time Sample Number')
    plt.ylabel('Power (Counts)')
    plt.title('Noise Diode Fold')
    if legend==True:
        plt.legend()


def plot_fullcalib(dio_cross,feedtype='l',**kwargs):
    """
    Generates and shows five plots: Uncalibrated diode, calibrated diode, fold information,
    phase offsets, and gain offsets for a noise diode measurement. Most useful diagnostic plot to
    make sure calibration proceeds correctly.
    """

    plt.figure("Multiple Calibration Plots", figsize=(12,9))
    left, width = 0.075,0.435
    bottom, height = 0.45,0.5
    width2 = 0.232
    bottom2, height2 = 0.115,0.0975

    rect_uncal = [left,bottom,width,height]
    rect_cal = [left+width+0.025,bottom,width,height]
    rect_fold = [left,bottom2,width2,0.22]
    rect_gain1 = [left+width2+0.1,bottom2,width2,height2]
    rect_phase1 = [left+width2*2+0.1*2,bottom2,width2,height2]
    rect_gain2 = [left+width2+0.1,bottom2+height2+0.025,width2,height2]
    rect_phase2 = [left+width2*2+0.1*2,bottom2+height2+0.025,width2,height2]

    #--------
    axFold = plt.axes(rect_fold)
    print('Plotting Diode Fold')
    plot_diode_fold(dio_cross,bothfeeds=False,feedtype=feedtype,min_samp=2000,max_samp=5500,legend=False,**kwargs)

    #--------
    print('Plotting Gain Offsets')
    plot_gain_offsets(dio_cross,feedtype=feedtype,ax1=rect_gain2,ax2=rect_gain1,legend=False,**kwargs)

    #--------
    print('Plotting Phase Offsets')
    plot_phase_offsets(dio_cross,feedtype=feedtype,ax1=rect_phase1,ax2=rect_phase2,legend=False,**kwargs)
    plt.ylabel('')

    #--------
    ax_uncal = plt.axes(rect_uncal)
    print('Plotting Uncalibrated Diode')
    plot_Stokes_diode(dio_cross,feedtype=feedtype,**kwargs)

    #--------
    ax_cal = plt.axes(rect_cal,sharey=ax_uncal)
    print('Plotting Calibrated Diode')
    plot_calibrated_diode(dio_cross,feedtype=feedtype,**kwargs)
    plt.ylabel('')
    plt.setp(ax_cal.get_yticklabels(),visible=False)

    plt.savefig(dio_cross[:-4]+'.stokescalib.png',dpi=2000)
    plt.show()

def plot_diodespec(ON_obs,OFF_obs,calflux,calfreq,spec_in,units='mJy',**kwargs):
    """
    Plots the full-band Stokes I spectrum of the noise diode (ON-OFF)
    """

    dspec = diode_spec(ON_obs,OFF_obs,calflux,calfreq,spec_in,**kwargs)
    obs = Waterfall(ON_obs,max_load=150)
    freqs = obs.populate_freqs()
    chan_per_coarse = obs.header['nchans']/obs.calc_n_coarse_chan()
    coarse_freqs = convert_to_coarse(freqs,chan_per_coarse)
    plt.ion()
    plt.figure()
    plt.plot(coarse_freqs,dspec)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Flux Density ('+units+')')
    plt.title('Noise Diode Spectrum')
