r"""
observatory.py
blimpy
"""

import os
import pandas as pd

DEBUGGING = False


class Observatory:
    r""" Class for handling observatory data.
        The Pandas dataframe is defined in ./observatory_info.csv
    """


    def __init__(self, telescope_id=None, telescope_name=None):
        """ init method for Observatory class
        Parameters (one or the other or neither):
            telescope_id (int): sigproc telescope_id
            telescope_name (str): telescope name
        If neither parameter is specified, then the fake telescope is used.
        Returns:
            self
        """
        abs_path = os.path.dirname(os.path.realpath(__file__))
        path_csv = os.path.join(abs_path, 'observatory_info.csv')
        df_full = pd.read_csv(path_csv, sep=',', engine='python', comment='#')
        if DEBUGGING:
            print("\nObservatory __init__ path_csv:", path_csv)
            print("Observatory __init__ df_full:\n", df_full)
        if telescope_id is not None:
            try:
                fields_dict = df_full.loc[df_full['TELESCOPE_ID'] == telescope_id].to_dict('list')
            except Exception as exc:
                raise RuntimeError("Observatory __init__ telescope_id={} is invalid!"
                                   .format(telescope_id)) from exc
        elif telescope_name is not None:
            try:
                fields_dict = df_full.loc[df_full['TELESCOPE_NAME'] == telescope_name].to_dict('list')
            except Exception as exc:
                raise RuntimeError("Observatory __init__ telescope_name={} is invalid!"
                                   .format(telescope_name)) from exc
        else:
            try:
                fields_dict = df_full.loc[df_full['TELESCOPE_ID'] == 0].to_dict('list')
            except Exception as exc:
                raise RuntimeError("Observatory __init__ Cannot find the fake telescope!") from exc
        if DEBUGGING:
            print("\nObservatory __init__ fields_dict:\n", fields_dict)
        self.telescope_name = fields_dict.get('TELESCOPE_NAME')[0]
        self.telescope_name_short = fields_dict.get('TELESCOPE_NAME_SHORT')[0]
        self.telescope_id = fields_dict.get('TELESCOPE_ID')[0]
        self.dish_diameter = fields_dict.get('DISH_DIAMETER')[0]
        self.xyz_coords = [
            fields_dict.get('X')[0],
            fields_dict.get('Y')[0],
            fields_dict.get('Z')[0] ]


    def get_telescope_name(self):
        r""" Return the telescope name to caller. """
        return self.telescope_name


    def get_telescope_name_short(self):
        r""" Return the short telescope name to caller. """
        return self.telescope_name_short


    def get_telescope_id(self):
        r""" Return the SIGPROC ID to caller. """
        return self.telescope_id


    def get_xyz_coords(self):
        r""" Return the X Y Z coordinates to caller. """
        return self.xyz_coords


    def get_dish_diameter(self):
        r""" Return the dish diameter to caller. """
        return self.dish_diameter


    def calc_beam_halfwidth(self, freq):
        """ Calculates beam halfwidth
            Code adapted from PRESTO
            Note: returns -1 if dish diameter = -1.
        Parameters:
            freq (int or float): frequency in MHz
        Returns:
            float: beam halfwidth in arcsec
        """
        # constants from PRESTO
        rad_to_deg = 57.29577951308232087679815481410517033240547246656
        sol = 299792458.0
        if self.dish_diameter == -1:
            return -1
        return 0.5 * 1.2 * sol / (freq * 1e6) / \
            self.dish_diameter * rad_to_deg * 3600.0


    def get_string(self):
        """ str method overload """
        fmt = 'Observatory: {}, short name: {}, telescope_id: {}, dish diameter: {}' \
            + ', (X, Y, Z): ({}, {}, {})'
        return fmt.format(self.telescope_name, self.telescope_name_short,
                          self.telescope_id, self.dish_diameter,
                          self.xyz_coords[0], self.xyz_coords[1], self.xyz_coords[2])


if __name__ == "__main__":
    obs = Observatory(telescope_id=6)
    print(obs.get_string())
