from blimpy import Waterfall
import numpy as np
from .fluxcal import foldcal

def get_stokes(cross_dat, feedtype='l'):
    '''Output stokes parameters (I,Q,U,V) for a rawspec
    cross polarization filterbank file'''

    #Compute Stokes Parameters
    if feedtype=='l':
        #I = XX+YY
        I = cross_dat[:,0,:]+cross_dat[:,1,:]
        #Q = XX-YY
        Q = cross_dat[:,0,:]-cross_dat[:,1,:]
        #U = 2*Re(XY)
        U = 2*cross_dat[:,2,:]
        #V = -2*Im(XY)
        V = -2*cross_dat[:,3,:]

    elif feedtype=='c':
        #I = LL+RR
        I = cross_dat[:,0,:]+cross_dat[:,1,:]
        #Q = 2*Re(RL)
        Q = 2*cross_dat[:,2,:]
        #U = 2*Im(RL)
        U = -2*cross_dat[:,3,:]
        #V = RR-LL
        V = cross_dat[:,1,:]-cross_dat[:,0,:]
    else:
        raise ValueError('feedtype must be \'l\' (linear) or \'c\' (circular)')

    #Add middle dimension to match Filterbank format
    I = np.expand_dims(I,axis=1)
    Q = np.expand_dims(Q,axis=1)
    U = np.expand_dims(U,axis=1)
    V = np.expand_dims(V,axis=1)

    #Compute linear polarization
    #L=np.sqrt(np.square(Q)+np.square(U))

    return I,Q,U,V

def convert_to_coarse(data,chan_per_coarse):
    '''
    Converts a data array with length n_chans to an array of length n_coarse_chans
    by averaging over the coarse channels
    '''
    #find number of coarse channels and reshape array
    num_coarse = data.size/chan_per_coarse
    data_shaped = np.array(np.reshape(data,(num_coarse,chan_per_coarse)))

    #Return the average over each coarse channel
    return np.mean(data_shaped[:,2:-1],axis=1)

def phase_offsets(Idat,Qdat,Udat,Vdat,tsamp,chan_per_coarse,feedtype='l',**kwargs):
    '''
    Calculates phase difference between X and Y feeds given U and V (U and Q for circular basis)
    data from a noise diode measurement on the target
    '''
    #Fold noise diode data and calculate ON OFF diferences for U and V
    if feedtype=='l':
        U_OFF,U_ON = foldcal(Udat,tsamp,**kwargs)
        V_OFF,V_ON = foldcal(Vdat,tsamp,**kwargs)
        Udiff = U_ON-U_OFF
        Vdiff = V_ON-V_OFF
        poffset = np.arctan2(-1*Vdiff,Udiff)

    if feedtype=='c':
        U_OFF,U_ON = foldcal(Udat,tsamp,**kwargs)
        Q_OFF,Q_ON = foldcal(Qdat,tsamp,**kwargs)
        Udiff = U_ON-U_OFF
        Qdiff = Q_ON-Q_OFF
        poffset = np.arctan2(Udiff,Qdiff)

    coarse_p =  convert_to_coarse(poffset,chan_per_coarse)

    #Correct for problems created by discontinuity in arctan
    #Find whether phase offsets have increasing or decreasing slope
    y = coarse_p[:6]
    x = np.arange(y.size)
    m = np.polyfit(x,y,1)[0]

    for i in range(coarse_p.size-3):
        if (m>0 and coarse_p[i+1]<coarse_p[i]) or (m<0 and coarse_p[i+1]>coarse_p[i]):
            coarse_p[i+1] = 2*coarse_p[i+2]-coarse_p[i+3]    #Move problem point near the next

    return coarse_p

def gain_offsets(Idat,Qdat,Udat,Vdat,tsamp,chan_per_coarse,feedtype='l',**kwargs):
    '''
    Determines relative gain error in the X and Y feeds for an
    observation given I and Q (I and V for circular basis) noise diode data.
    '''
    if feedtype=='l':
        #Fold noise diode data and calculate ON OFF differences for I and Q
        I_OFF,I_ON = foldcal(Idat,tsamp,**kwargs)
        Q_OFF,Q_ON = foldcal(Qdat,tsamp,**kwargs)

        #Calculate power in each feed for noise diode ON and OFF
        XX_OFF = (I_OFF+Q_OFF)/2
        YY_OFF = (I_OFF-Q_OFF)/2

        #Calculate gain offset (divided by 2) as defined in Heiles (2001)
        G = (XX_OFF-YY_OFF)/(XX_OFF+YY_OFF)

    if feedtype=='c':
        #Fold noise diode data and calculate ON OFF differences for I and Q
        I_OFF,I_ON = foldcal(Idat,tsamp,**kwargs)
        V_OFF,V_ON = foldcal(Vdat,tsamp,**kwargs)

        #Calculate power in each feed for noise diode ON and OFF
        RR_OFF = (I_OFF+V_OFF)/2
        LL_OFF = (I_OFF-V_OFF)/2

        #Calculate gain offset (divided by 2) as defined in Heiles (2001)
        G = (RR_OFF-LL_OFF)/(RR_OFF+LL_OFF)

    return convert_to_coarse(G,chan_per_coarse)

def apply_Mueller(I,Q,U,V, gain_offsets, phase_offsets, chan_per_coarse, feedtype='l'):
    '''
    Returns calibrated Stokes parameters for an observation given an array
    of differential gains and phase differences.
    '''

    #Find shape of data arrays and calculate number of coarse channels
    shape = I.shape
    ax0 = I.shape[0]
    ax1 = I.shape[1]
    nchans = I.shape[2]
    ncoarse = nchans/chan_per_coarse

    #Reshape data arrays to separate coarse channels
    I = np.reshape(I,(ax0,ax1,ncoarse,chan_per_coarse))
    Q = np.reshape(Q,(ax0,ax1,ncoarse,chan_per_coarse))
    U = np.reshape(U,(ax0,ax1,ncoarse,chan_per_coarse))
    V = np.reshape(V,(ax0,ax1,ncoarse,chan_per_coarse))

    #Swap axes 2 and 3 to in order for broadcasting to work correctly
    I = np.swapaxes(I,2,3)
    Q = np.swapaxes(Q,2,3)
    U = np.swapaxes(U,2,3)
    V = np.swapaxes(V,2,3)

    #Apply top left corner of electronics chain inverse Mueller matrix
    a = 1/(1-gain_offsets**2)
    if feedtype=='l':
        Icorr = a*(I-gain_offsets*Q)
        Qcorr = a*(-1*gain_offsets*I+Q)
        I = None
        Q = None
    if feedtype=='c':
        Icorr = a*(I-gain_offsets*V)
        Vcorr = a*(-1*gain_offsets*I+V)
        I = None
        V = None

    #Apply bottom right corner of electronics chain inverse Mueller matrix
    if feedtype=='l':
        Ucorr = U*np.cos(phase_offsets)-V*np.sin(phase_offsets)
        Vcorr = U*np.sin(phase_offsets)+V*np.cos(phase_offsets)
        U = None
        V = None
    if feedtype=='c':
        Qcorr = Q*np.cos(phase_offsets)+U*np.sin(phase_offsets)
        Ucorr = -1*Q*np.sin(phase_offsets)+U*np.cos(phase_offsets)
        Q = None
        U = None

    #Reshape arrays to original shape
    Icorr = np.reshape(np.swapaxes(Icorr,2,3),shape)
    Qcorr = np.reshape(np.swapaxes(Qcorr,2,3),shape)
    Ucorr = np.reshape(np.swapaxes(Ucorr,2,3),shape)
    Vcorr = np.reshape(np.swapaxes(Vcorr,2,3),shape)

    #Return corrected data arrays
    return Icorr,Qcorr,Ucorr,Vcorr

def calibrate_pols(cross_pols,diode_cross,obsI=None,onefile=True,feedtype='l',**kwargs):
    '''
    Write Stokes-calibrated filterbank file for a given observation
    with a calibrator noise diode measurement on the source

    Parameters
    ----------
    cross_pols : string
        Path to cross polarization filterbank file (rawspec output) for observation to be calibrated
    diode_cross : string
        Path to cross polarization filterbank file of noise diode measurement ON the target
    obsI : string
        Path to Stokes I filterbank file of main observation (only needed if onefile=False)
    onefile : boolean
        True writes all calibrated Stokes parameters to a single filterbank file,
        False writes four separate files
    feedtype : 'l' or 'c'
        Basis of antenna dipoles. 'c' for circular, 'l' for linear
    '''
    #Obtain time sample length, frequencies, and noise diode data
    obs = Waterfall(diode_cross,max_load=150)
    cross_dat = obs.data
    tsamp = obs.header['tsamp']

    #Calculate number of coarse channels in the noise diode measurement (usually 8)
    dio_ncoarse = obs.calc_n_coarse_chan()
    dio_nchans = obs.header['nchans']
    dio_chan_per_coarse = dio_nchans/dio_ncoarse
    obs = None
    Idat,Qdat,Udat,Vdat = get_stokes(cross_dat,feedtype)
    cross_dat = None
    #Calculate differential gain and phase from noise diode measurements
    print('Calculating Mueller Matrix variables')
    gams = gain_offsets(Idat,Qdat,Udat,Vdat,tsamp,dio_chan_per_coarse,feedtype,**kwargs)
    psis = phase_offsets(Idat,Qdat,Udat,Vdat,tsamp,dio_chan_per_coarse,feedtype,**kwargs)

    #Clear data arrays to save memory
    Idat = None
    Qdat = None
    Udat = None
    Vdat = None

    #Get corrected Stokes parameters
    print('Opening '+cross_pols)
    cross_obs = Waterfall(cross_pols,max_load=150)
    obs_ncoarse = cross_obs.calc_n_coarse_chan()
    obs_nchans = cross_obs.header['nchans']
    obs_chan_per_coarse = obs_nchans/obs_ncoarse

    print('Grabbing Stokes parameters')
    I,Q,U,V = get_stokes(cross_obs.data,feedtype)

    print('Applying Mueller Matrix')
    I,Q,U,V = apply_Mueller(I,Q,U,V,gams,psis,obs_chan_per_coarse,feedtype)

    #Use onefile (default) to produce one filterbank file containing all Stokes information
    if onefile:
        cross_obs.data[:,0,:] = np.squeeze(I)
        cross_obs.data[:,1,:] = np.squeeze(Q)
        cross_obs.data[:,2,:] = np.squeeze(U)
        cross_obs.data[:,3,:] = np.squeeze(V)
        cross_obs.write_to_fil(cross_pols[:-15]+'.SIQUV.polcal.fil')
        print('Calibrated Stokes parameters written to '+cross_pols[:-15]+'.SIQUV.polcal.fil')
        return

    # If onefile=False and obsI=None, abort
    assert obsI is not None

    #Write corrected Stokes parameters to four filterbank files if onefile==False
    obs = Waterfall(obsI,max_load=150)
    obs.data = I
    obs.write_to_fil(cross_pols[:-15]+'.SI.polcal.fil')   #assuming file is named *.cross_pols.fil
    print('Calibrated Stokes I written to '+cross_pols[:-15]+'.SI.polcal.fil')

    obs.data = Q
    obs.write_to_fil(cross_pols[:-15]+'.Q.polcal.fil')   #assuming file is named *.cross_pols.fil
    print('Calibrated Stokes Q written to '+cross_pols[:-15]+'.Q.polcal.fil')

    obs.data = U
    obs.write_to_fil(cross_pols[:-15]+'.U.polcal.fil')   #assuming file is named *.cross_pols.fil
    print('Calibrated Stokes U written to '+cross_pols[:-15]+'.U.polcal.fil')

    obs.data = V
    obs.write_to_fil(cross_pols[:-15]+'.V.polcal.fil')   #assuming file is named *.cross_pols.fil
    print('Calibrated Stokes V written to '+cross_pols[:-15]+'.V.polcal.fil')


def fracpols(cross_dat, **kwargs):
    '''Output fractional linear and circular polarizations for a
    rawspec cross polarization .fil file. NOT STANDARD USE'''

    I,Q,U,V,L=get_stokes(cross_dat, **kwargs)
    return L/I,V/I

def write_stokefils(cross_dat, str_I, Ifil=False, Qfil=False, Ufil=False, Vfil=False, Lfil=False, **kwargs):
    '''Writes up to 5 new filterbank files corresponding to each Stokes
    parameter (and total linear polarization L) for a given cross polarization .fil file'''

    I,Q,U,V,L=get_stokes(cross_dat, **kwargs)
    obs = Waterfall(str_I, max_load=150) #Load filterbank file to write stokes data to
    if Ifil:
        obs.data = I
        obs.write_to_fil(cross_dat[:-15]+'.I.fil')   #assuming file is named *.cross_pols.fil

    if Qfil:
        obs.data = Q
        obs.write_to_fil(cross_dat[:-15]+'.Q.fil')   #assuming file is named *.cross_pols.fil

    if Ufil:
        obs.data = U
        obs.write_to_fil(cross_dat[:-15]+'.U.fil')   #assuming file is named *.cross_pols.fil

    if Vfil:
        obs.data = V
        obs.write_to_fil(cross_dat[:-15]+'.V.fil')   #assuming file is named *.cross_pols.fil

    if Lfil:
        obs.data = L
        obs.write_to_fil(cross_dat[:-15]+'.L.fil')   #assuming file is named *.cross_pols.fil


def write_polfils(cross_dat, str_I, **kwargs):
    '''Writes two new filterbank files containing fractional linear and
    circular polarization data'''

    lin,circ=fracpols(cross_dat, **kwargs)
    obs = Waterfall(str_I, max_load=150)

    obs.data = lin
    obs.write_to_fil(cross_dat[:-15]+'.linpol.fil')   #assuming file is named *.cross_pols.fil

    obs.data = circ
    obs.write_to_fil(cross_dat[:-15]+'.circpol.fil')   #assuming file is named *.cross_pols.fil
