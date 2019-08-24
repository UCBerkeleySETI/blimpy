from blimpy import Waterfall
import numpy as np
import pylab as plt
from scipy.integrate import trapz
from astropy import units as u
import matplotlib.pyplot as plt

def foldcal(data,tsamp, diode_p=0.04,numsamps=1000,switch=False,inds=False):
    '''
    Returns time-averaged spectra of the ON and OFF measurements in a
    calibrator measurement with flickering noise diode

    Parameters
    ----------
    data : 2D Array object (float)
        2D dynamic spectrum for data (any Stokes parameter) with flickering noise diode.
    tsamp : float
        Sampling time of data in seconds
    diode_p : float
        Period of the flickering noise diode in seconds
    numsamps : int
        Number of samples over which to average noise diode ON and OFF
    switch : boolean
        Use switch=True if the noise diode "skips" turning from OFF to ON once or vice versa
    inds : boolean
        Use inds=True to also return the indexes of the time series where the ND is ON and OFF
    '''

    halfper = diode_p/2.0

    foldt = halfper/tsamp   #number of time samples per diode switch

    onesec = 1/tsamp    #number of time samples in the first second

    #Find diode switches in units of time samples and round down to the nearest int
    ints = np.arange(0,numsamps)
    t_switch = (onesec+ints*foldt)
    t_switch = t_switch.astype('int')

    ONints = np.array(np.reshape(t_switch[:],(numsamps/2,2)))
    ONints[:,0] = ONints[:,0]+1   #Find index ranges of ON time samples

    OFFints = np.array(np.reshape(t_switch[1:-1],(numsamps/2-1,2)))
    OFFints[:,0] = OFFints[:,0]+1   #Find index ranges of OFF time samples

    av_ON = []
    av_OFF = []

    #Average ON and OFF spectra separately with respect to time
    for i in ONints:
        if i[1]!=i[0]:
            av_ON.append(np.sum(data[i[0]:i[1],:,:],axis=0)/(i[1]-i[0]))

    for i in OFFints:
        if i[1]!=i[0]:
            av_OFF.append(np.sum(data[i[0]:i[1],:,:],axis=0)/(i[1]-i[0]))

    #If switch=True, flip the return statement since ON is actually OFF
    if switch==False:
        if inds==False:
            return np.squeeze(np.mean(av_ON,axis=0)), np.squeeze(np.mean(av_OFF,axis=0))
        else:
            return np.squeeze(np.mean(av_ON,axis=0)), np.squeeze(np.mean(av_OFF,axis=0)),ONints,OFFints
    if switch==True:
        if inds==False:
            return np.squeeze(np.mean(av_OFF,axis=0)), np.squeeze(np.mean(av_ON,axis=0))
        else:
            return np.squeeze(np.mean(av_OFF,axis=0)), np.squeeze(np.mean(av_ON,axis=0)),OFFints,ONints

def integrate_chans(spec,freqs,chan_per_coarse):
    '''
    Integrates over each core channel of a given spectrum.
    Important for calibrating data with frequency/time resolution different from noise diode data

    Parameters
    ----------
    spec : 1D Array (float)
        Spectrum (any Stokes parameter) to be integrated
    freqs : 1D Array (float)
        Frequency values for each bin of the spectrum
    chan_per_coarse: int
        Number of frequency bins per coarse channel
    '''

    num_coarse = spec.size/chan_per_coarse     #Calculate total number of coarse channels

    #Rearrange spectrum by coarse channel
    spec_shaped = np.array(np.reshape(spec,(num_coarse,chan_per_coarse)))
    freqs_shaped = np.array(np.reshape(freqs,(num_coarse,chan_per_coarse)))

    #Average over coarse channels
    return np.mean(spec_shaped[:,1:-1],axis=1)

def integrate_calib(name,chan_per_coarse,fullstokes=False,**kwargs):
    '''
    Folds Stokes I noise diode data and integrates along coarse channels

    Parameters
    ----------
    name : str
        Path to noise diode filterbank file
    chan_per_coarse : int
        Number of frequency bins per coarse channel
    fullstokes : boolean
        Use fullstokes=True if data is in IQUV format or just Stokes I, use fullstokes=False if
        it is in cross_pols format
    '''
    #Load data
    obs = Waterfall(name,max_load=150)
    data = obs.data

    #If the data has cross_pols format calculate Stokes I
    if fullstokes==False and data.shape[1]>1:
        data = data[:,0,:]+data[:,1,:]
        data = np.expand_dims(data,axis=1)
    #If the data has IQUV format get Stokes I
    if fullstokes==True:
        data = data[:,0,:]
        data = np.expand_dims(data,axis=1)

    tsamp = obs.header['tsamp']

    #Calculate ON and OFF values
    OFF,ON = foldcal(data,tsamp,**kwargs)

    freqs = obs.populate_freqs()

    #Find ON and OFF spectra by coarse channel
    ON_int = integrate_chans(ON,freqs,chan_per_coarse)
    OFF_int = integrate_chans(OFF,freqs,chan_per_coarse)

    #If "ON" is actually "OFF" switch them
    if np.sum(ON_int)<np.sum(OFF_int):
        temp = ON_int
        ON_int = OFF_int
        OFF_int = temp

    #Return coarse channel spectrum of OFF and ON
    return OFF_int,ON_int

def get_calfluxes(calflux,calfreq,spec_in,centerfreqs,oneflux):
    '''
    Given properties of the calibrator source, calculate fluxes of the source
    in a particular frequency range

    Parameters
    ----------
    calflux : float
        Known flux of calibrator source at a particular frequency
    calfreq : float
        Frequency where calibrator source has flux calflux (MHz) (see above)
    spec_in : float
        Known power-law spectral index of calibrator source. Use convention flux(frequency) = constant * frequency^(spec_in)
    centerfreqs : 1D Array (float)
        Central frequency values of each coarse channel
    oneflux : boolean
        Use oneflux to choose between calculating the flux for each core channel (False)
        or using one value for the entire frequency range (True)
    '''

    const = calflux/np.power(calfreq,spec_in)
    if oneflux==False:
        return const*np.power(centerfreqs,spec_in)
    else:
        return const*np.power(np.mean(centerfreqs),spec_in)

def get_centerfreqs(freqs,chan_per_coarse):
    '''
    Returns central frequency of each coarse channel

    Parameters
    ----------
    freqs : 1D Array (float)
        Frequency values for each bin of the spectrum
    chan_per_coarse: int
        Number of frequency bins per coarse channel
    '''

    num_coarse = freqs.size/chan_per_coarse
    freqs = np.reshape(freqs,(num_coarse,chan_per_coarse))
    return np.mean(freqs,axis=1)

def f_ratios(calON_obs,calOFF_obs,chan_per_coarse,**kwargs):
    '''
    Calculate f_ON, and f_OFF as defined in van Straten et al. 2012 equations 2 and 3

    Parameters
    ----------
    calON_obs : str
        Path to filterbank file (any format) for observation ON the calibrator source
    calOFF_obs : str
        Path to filterbank file (any format) for observation OFF the calibrator source
    '''
    #Calculate noise diode ON and noise diode OFF spectra (H and L) for both observations
    L_ON,H_ON = integrate_calib(calON_obs,chan_per_coarse,**kwargs)
    L_OFF,H_OFF = integrate_calib(calOFF_obs,chan_per_coarse,**kwargs)

    f_ON = H_ON/L_ON-1
    f_OFF = H_OFF/L_OFF-1

    return f_ON, f_OFF


def diode_spec(calON_obs,calOFF_obs,calflux,calfreq,spec_in,average=True,oneflux=False,**kwargs):
    '''
    Calculate the coarse channel spectrum and system temperature of the noise diode in Jy given two noise diode
    measurements ON and OFF the calibrator source with the same frequency and time resolution

    Parameters
    ----------
    calON_obs : str
        (see f_ratios() above)
    calOFF_obs : str
        (see f_ratios() above)
    calflux : float
        Known flux of calibrator source at a particular frequency
    calfreq : float
        Frequency where calibrator source has flux calflux (see above)
    spec_in : float
        Known power-law spectral index of calibrator source. Use convention flux(frequency) = constant * frequency^(spec_in)
    average : boolean
        Use average=True to return noise diode and Tsys spectra averaged over frequencies
    '''
    #Load frequencies and calculate number of channels per coarse channel
    obs = Waterfall(calON_obs,max_load=150)
    freqs = obs.populate_freqs()
    ncoarse = obs.calc_n_coarse_chan()
    nchans = obs.header['nchans']
    chan_per_coarse = nchans/ncoarse

    f_ON, f_OFF = f_ratios(calON_obs,calOFF_obs,chan_per_coarse,**kwargs)

    #Obtain spectrum of the calibrator source for the given frequency range
    centerfreqs = get_centerfreqs(freqs,chan_per_coarse)
    calfluxes = get_calfluxes(calflux,calfreq,spec_in,centerfreqs,oneflux)

    #C_o and Tsys as defined in van Straten et al. 2012
    C_o = calfluxes/(1/f_ON-1/f_OFF)
    Tsys = C_o/f_OFF

    #return coarse channel diode spectrum
    if average==True:
        return np.mean(C_o),np.mean(Tsys)
    else:
        return C_o,Tsys

def get_Tsys(calON_obs,calOFF_obs,calflux,calfreq,spec_in,oneflux=False,**kwargs):
    '''
    Returns frequency dependent system temperature given observations on and off a calibrator source

    Parameters
    ----------
    (See diode_spec())
    '''
    obs = Waterfall(calON_obs,max_load=150)
    ncoarse = obs.calc_n_coarse_chan()
    chan_per_coarse = obs.header['nchans']/ncoarse
    freqs = obs.populate_freqs()
    cfreqs = get_centerfreqs(freqs,chan_per_coarse)
    S_sys = diode_spec(calON_obs,calOFF_obs,calflux,calfreq,spec_in,average=False,oneflux=False,**kwargs)[1]

    T_sys = Jy_to_Kelvin(S_sys,cfreqs)
    return T_sys,cfreqs

def get_Tsys_nodiode(calON_obs_name,calOFF_obs_name,calflux,calfreq,spec_in,plot=False,**kwargs):
    '''
    Calculates system temperature from two flux calibrator scans taken without noise diode flickering.
    CURRENTLY ONLY IMPLEMENTED FOR STOKES I DATA.

    Parameters
    ----------
    calON_obs_name : str
        Path to filterbank file for scan ON calibrator target
    calOFF_obs_name : str
        Path to filterbank file for scan OFF calibrator
    calflux : float
        Flux in Jy of the calibrator source
    calfreq : float or int
        Frequency at which calflux was taken (MHz)
    spec_in : float
        Spectral index of this calibrator
    '''
    calON_obs = Waterfall(calON_obs_name,max_load=150)
    calOFF_obs = Waterfall(calOFF_obs_name,max_load=150)
    ncoarse = calON_obs.calc_n_coarse_chan()
    chan_per_coarse = calON_obs.header['nchans']/ncoarse

    ONobs_spec = np.squeeze(np.mean(calON_obs.data,axis=0))
    OFFobs_spec = np.squeeze(np.mean(calOFF_obs.data,axis=0))
    freqs = calON_obs.populate_freqs()

    ON_coarse_spec = integrate_chans(ONobs_spec,freqs,chan_per_coarse)
    OFF_coarse_spec = integrate_chans(OFFobs_spec,freqs,chan_per_coarse)
    cfreqs = get_centerfreqs(freqs,chan_per_coarse)
    cal_fluxes = get_calfluxes(calflux,calfreq,spec_in,cfreqs,oneflux=False)

    S_sys = (OFF_coarse_spec)/(ON_coarse_spec-OFF_coarse_spec)*cal_fluxes

    #Convert to Kelvins
    T_sys,cfreqs = Jy_to_Kelvin(S_sys,cfreqs)
    return T_sys,cfreqs


def Jy_to_Kelvin(flux,freqs):
    '''
    Convert flux spectrum in Jy/beam to temperature units. Frequency inputs in MHz.

    Parameters
    ----------
    flux : 1D Array (float)
        Spectrum over frequency band in Jy
    freqs : 1D Array (float)
        Frequencies (in MHz)
    '''
    beam_FWHM = 755.0/(np.mean(freqs*10**(-3)))*u.arcsec
    fwhm_to_sig = 1./(8*np.log(2))**0.5
    beam_area = 2*np.pi*(beam_FWHM*fwhm_to_sig)**2
    freq_MHz = np.mean(freqs)*u.MHz
    equiv = u.brightness_temperature(beam_area,freq_MHz)
    T = (flux*u.Jy).to(u.K,equivalencies=equiv)

    return T,freqs

def calibrate_fluxes(main_obs_name,dio_name,dspec,Tsys,fullstokes=False,**kwargs):
    '''
    Produce calibrated Stokes I for an observation given a noise diode
    measurement on the source and a diode spectrum with the same number of
    coarse channels

    Parameters
    ----------
    main_obs_name : str
        Path to filterbank file containing final data to be calibrated
    dio_name : str
        Path to filterbank file for observation on the target source with flickering noise diode
    dspec : 1D Array (float) or float
        Coarse channel spectrum (or average) of the noise diode in Jy (obtained from diode_spec())
    Tsys : 1D Array (float) or float
        Coarse channel spectrum (or average) of the system temperature in Jy
    fullstokes: boolean
        Use fullstokes=True if data is in IQUV format or just Stokes I, use fullstokes=False if
        it is in cross_pols format
    '''

    #Find folded spectra of the target source with the noise diode ON and OFF
    main_obs = Waterfall(main_obs_name,max_load=150)
    ncoarse = main_obs.calc_n_coarse_chan()
    dio_obs = Waterfall(dio_name,max_load=150)
    dio_chan_per_coarse = dio_obs.header['nchans']/ncoarse
    dOFF,dON = integrate_calib(dio_name,dio_chan_per_coarse,fullstokes,**kwargs)

    #Find Jy/count for each coarse channel using the diode spectrum
    main_dat = main_obs.data
    scale_facs = dspec/(dON-dOFF)
    print(scale_facs)

    nchans = main_obs.header['nchans']
    obs_chan_per_coarse = nchans/ncoarse

    ax0_size = np.size(main_dat,0)
    ax1_size = np.size(main_dat,1)

    #Reshape data array of target observation and multiply coarse channels by the scale factors
    main_dat = np.reshape(main_dat,(ax0_size,ax1_size,ncoarse,obs_chan_per_coarse))
    main_dat = np.swapaxes(main_dat,2,3)

    main_dat = main_dat*scale_facs
    main_dat = main_dat-Tsys
    main_dat = np.swapaxes(main_dat,2,3)
    main_dat = np.reshape(main_dat,(ax0_size,ax1_size,nchans))

    #Write calibrated data to a new filterbank file with ".fluxcal" extension

    main_obs.data = main_dat
    main_obs.write_to_filterbank(main_obs_name[:-4]+'.fluxcal.fil')
    print('Finished: calibrated product written to ' + main_obs_name[:-4]+'.fluxcal.fil')


#end module
