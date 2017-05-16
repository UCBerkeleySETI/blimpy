#!/usr/bin/env python
"""
setup.py -- setup script for use of packages.
"""
from setuptools import setup, find_packages

version = '1.1.1'

# create entry points
# see http://astropy.readthedocs.org/en/latest/development/scripts.html
entry_points = {
    'console_scripts' :
        ['filutil = blimpy.filterbank:cmd_tool',
         'watutil = blimpy.waterfall:cmd_tool',
         'rawutil = blimpy.guppi:cmd_tool',
         'fil2hdf = blimpy.fil2hdf:cmd_tool',
         'gup2hdf = blimpy.gup2hdf:cmd_tool',
         'fil2h5 = blimpy.fil2h5:make_h5_file'
     ]
    }

setup(name='blimpy',
      version=version,
      description='Python utilities for Breakthrough Listen SETI observations',
      install_requires=['astropy', 'numpy', 'cython', 'h5py'],
      url='https://github.com/ucberkeleyseti/blimpy',
      author='Danny Price',
      author_email='dancpr@berkeley.edu',
      license='MIT',
      entry_points=entry_points,
      packages=find_packages(),
      zip_safe=False,
      )
