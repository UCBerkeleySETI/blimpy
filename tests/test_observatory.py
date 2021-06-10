r"""
test_observatory.py
"""

from blimpy.ephemeris import Observatory


def error_msg(arg_string):
    r""" Just making clearer error messages """
    return "test_observatory.py: " + arg_string


def test_observatory_construction():
    r""" Constructor test """
    print("\n===== Begin test_observatory_construction")
    obs = Observatory()
    assert obs.get_telescope_name() == "Fake", error_msg("Wrong name for the fake observatory")
    obs = Observatory(telescope_id=4)
    assert obs.get_telescope_name() == "PARKES", error_msg("Wrong name for the Parkes observatory")
    assert obs.get_telescope_name_short() == "PK", \
        error_msg("Wrong short name for the Parkes observatory")
    obs = Observatory(telescope_name="GBT")
    assert obs.get_telescope_id() == 6, error_msg("Wrong Sigproc ID for the GBT observatory")
    print("===== End test_observatory_construction")


def test_observatory_values():
    r""" Observatory values test along with beam halfwidth calculation test"""
    print("\n===== Begin test_observatory_values")
    obs = Observatory(telescope_id=0)
    print(obs.get_string())
    assert obs.get_telescope_name() == 'Fake', error_msg("Incorrect telescope name")
    assert obs.get_xyz_coords() == [0.0, 0.0, 0.0], error_msg("Incorrect XYZ coords")

    gbt = Observatory(telescope_id=6)
    beam_halfwidth = gbt.calc_beam_halfwidth(100)
    assert (beam_halfwidth - 3710.19799582) < .0000001, \
        error_msg("Incorrect beam haflwidth calculation")
    print("===== End test_observatory_values")


def test_observatory_procs():
    r""" Try the member functions with Parkes """
    print("\n===== Begin test_observatory_procs")
    obs = Observatory(telescope_id=4)
    print(obs.get_string())
    print(obs.get_telescope_name())
    print(obs.get_telescope_name_short())
    print(obs.get_telescope_id())
    print(obs.get_xyz_coords())
    print(obs.get_dish_diameter())
    print(obs.get_string())
    print("===== End test_observatory_procs")


if __name__ == "__main__":
    test_observatory_construction()
    test_observatory_values()
    test_observatory_procs()
