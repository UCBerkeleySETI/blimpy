#!/usr/local/bin/python
"""-------------------------------------
observatory.py
blimpy
--------------------------------------
"""
import os
import sys
import pandas as pd


class Observatory:
    """ Class for handling observatory data.
        Data is stored in ./observatory_info.csv
    """

    def __init__(self, telescope_id=None, telescope_name=None):
        """ init method for Observatory class
        Parameters:
            telescope_id (int): sigproc telescope_id
            telescope_name (str): telescope name. see ./observatory_info.csv
        Returns:
            self
        """
        abs_path = os.path.dirname(os.path.realpath(__file__))
        iter_csv = pd.read_csv(os.path.join(abs_path, 'observatory_info.csv'),
            comment='#',
            iterator=True,
            chunksize=1000)
      
        if telescope_id is not None:
            dataframe = pd.concat(
                [chunk[chunk['SIGPROC_ID'] == telescope_id] for chunk in iter_csv])
        elif telescope_name is not None:
            dataframe = pd.concat(
                [chunk[chunk['TELESCOPE_NAME'] == telescope_name] for chunk in iter_csv])
        else:
            dataframe = pd.concat(
                [chunk[chunk['SIGPROC_ID'] == 0] for chunk in iter_csv])
        fields_dict = dataframe.to_dict("record")[0]
        self._telescope_name = fields_dict['TELESCOPE_NAME']
        self._telescope_name_short = fields_dict['TELESCOPE_NAME_SHORT']
        self._sigproc_id = fields_dict['SIGPROC_ID']
        self._xyz_coords = (
            float(fields_dict['X']),
            float(fields_dict['Y']),
            float(fields_dict['Z']))
        self._dish_diam = fields_dict['DISH_DIAMETER']
        self._info_dict = fields_dict

    def get_telescope_name(self):
        """getter method"""
        return self._telescope_name

    def get_telescope_name_short(self):
        """getter method"""
        return self._telescope_name_short

    def get_sigproc_id(self):
        """getter method"""
        return self._sigproc_id

    def get_xyz_coords(self):
        """getter method"""
        return self._xyz_coords

    def calc_beam_halfwidth(self, freq):
        """ Calculates beam halfwidth
            Code adapted from PRESTO
            Note: returns -1 if dish diameter. Possible //TODO, throw error
        Parameters:
            freq (int or float): frequency in MHz
        Returns:
            float: beam halfwidth in arcsec
        """
        # constants from PRESTO
        rad_to_deg = 57.29577951308232087679815481410517033240547246656
        sol = 299792458.0
        if self._dish_diam == -1:
            return -1
        return 0.5 * 1.2 * sol / (freq * 1e6) / \
            self._dish_diam * rad_to_deg * 3600.0

    def __str__(self):
        """ str method overload """
        output_str = 'Observatory: ' + self._telescope_name + "\n"
        keys = sorted(self._info_dict.keys())
        for key in keys:
            if key != "TELESCOPE_NAME" and "ADDED" not in key:
                output_str += "\t" + str(key) + ": " + \
                    str(self._info_dict[key]) + "\n"
        return output_str
