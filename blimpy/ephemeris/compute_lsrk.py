from .config import *
from .compute_lst import compute_lst

def compute_lsrk(wf):
    """ Computes the LSR in km/s

    Computes the Local standard of rest kinematic using the time (MJD),
    RA and DEC of the observation to compute along with the telescope location.
    Requires pyslalib

    Args:
        wf (bl.Waterfall): Waterfall object for which to compute LSR
    """
    ra = Angle(wf.header['src_raj'], unit='hourangle')
    dec = Angle(wf.header['src_dej'], unit='degree')
    mjdd = wf.header['tstart']
    rarad = ra.to('radian').value
    dcrad = dec.to('radian').value
    last = compute_lst(wf)
    tellat = np.deg2rad(wf.coords[0])
    tellong = np.deg2rad(wf.coords[1])

    # convert star position to vector
    starvect = s.sla_dcs2c(rarad, dcrad)

    # velocity component in ra,dec due to Earth rotation
    Rgeo = s.sla_rverot(tellat, rarad, dcrad, last)

    # get Barycentric and heliocentric velocity and position of the Earth.
    evp = s.sla_evp(mjdd, 2000.0)
    dvb = evp[0]  # barycentric velocity vector, in AU/sec
    dpb = evp[1]  # barycentric position vector, in AU
    dvh = evp[2]  # heliocentric velocity vector, in AU/sec
    dph = evp[3]  # heliocentric position vector, in AU

    # dot product of vector to object and heliocentric velocity
    # convert AU/sec to km/sec
    vcorhelio = -s.sla_dvdv(starvect, dvh) * 149.597870e6
    vcorbary = -s.sla_dvdv(starvect, dvb) * 149.597870e6

    # rvlsrd is velocity component in ra,dec direction due to the Sun's
    # motion with respect to the "dynamical" local standard of rest
    rvlsrd = s.sla_rvlsrd(rarad, dcrad)

    # rvlsrk is velocity component in ra,dec direction due to i
    # the Sun's motion w.r.t the "kinematic" local standard of rest
    rvlsrk = s.sla_rvlsrk(rarad, dcrad)

    # rvgalc is velocity component in ra,dec direction due to
    # the rotation of the Galaxy.
    rvgalc = s.sla_rvgalc(rarad, dcrad)
    totalhelio = Rgeo + vcorhelio
    totalbary = Rgeo + vcorbary
    totallsrk = totalhelio + rvlsrk
    totalgal = totalbary + rvlsrd + rvgalc

    return totallsrk
