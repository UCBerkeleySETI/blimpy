from .config import *

def compute_lst(wf):
    """ Compute LST for observation

    Computes local sidereal time (LST) for the observation, using SLALIB.

    Args:
        wf (bl.Waterfall): blimpy Waterfall object.
    """
    if wf.header['telescope_id'] == 6:
        wf.coords = gbt_coords
    elif wf.header['telescope_id'] == 4:
        wf.coords = parkes_coords
    else:
        raise RuntimeError("Currently only Parkes and GBT supported")
    if HAS_SLALIB:
        # dut1 = (0.2 /3600.0) * np.pi/12.0
        dut1 = 0.0
        mjd = wf.header['tstart']
        tellong = np.deg2rad(wf.coords[1])
        last = s.sla_gmst(mjd) - tellong + s.sla_eqeqx(mjd) + dut1
        # lmst = s.sla_gmst(mjd) - tellong
        if last < 0.0: last = last + 2.0 * np.pi
        return last
    else:
        raise RuntimeError("This method requires pySLALIB")
