"""
test_heavy.py
BlimPy

"""

from blimpy.ephemeris import Observatory

def error_msg(s):
    """ Just making clearer error messages """
    return "test_observatory.py: " + s

def test_observatory_construction():
    """ Constructor test """
    obs = Observatory(telescope_id=0)
    assert obs.get_telescope_name() != None, error_msg("Could not create observatory")

def test_observatory_values():
    """ Observatory values test along with beam halfwidth calculation test"""
    obs = Observatory(telescope_id=0)
  
    assert obs.get_telescope_name() == 'Fake', error_msg("Incorrect name")
    assert obs.get_xyz_coords() == (0.0,0.0,0.0), error_msg("Incorrect XYZ coords")

    gbt = Observatory(telescope_id=6)
    beam_halfwidth = gbt.calc_beam_halfwidth(100)

    assert (beam_halfwidth - 3710.19799582) < .0000001, error_msg("Incorrect beam haflwidth calculation")

if __name__ == "__main__":
    test_observatory_construction()
    test_observatory_values()
