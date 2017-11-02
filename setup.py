#!/usr/bin/env python
"""
setup.py -- setup script for use of packages.
"""
from setuptools import setup, find_packages

version = '1.1.6'

# create entry points
# see http://astropy.readthedocs.org/en/latest/development/scripts.html
entry_points = {
    'console_scripts' :
        ['filutil = blimpy.filterbank:cmd_tool',
         'watutil = blimpy.waterfall:cmd_tool',
         'rawutil = blimpy.guppi:cmd_tool',
#         'fil2hdf = blimpy.fil2hdf:cmd_tool',  #EE deprecating until tested
#         'gup2hdf = blimpy.gup2hdf:cmd_tool',  #EE deprecating until tested
         'fil2h5 = blimpy.fil2h5:cmd_tool',
         'h52fil = blimpy.h52fil:cmd_tool',
         'matchfils = blimpy.match_fils:cmd_tool'
     ]
    }

setup(name='blimpy',
      version=version,
      description='Python utilities for Breakthrough Listen SETI observations',
      install_requires=['astropy', 'numpy', 'cython', 'h5py'],
      url='https://github.com/ucberkeleyseti/blimpy',
      author='Danny Price & Emilio Enriquez',
      author_email='dancpr@berkeley.edu',
      license='MIT',
      entry_points=entry_points,
      packages=find_packages(),
      zip_safe=False,
      )
