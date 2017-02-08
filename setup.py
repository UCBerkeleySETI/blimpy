#!/usr/bin/env python
"""
setup.py -- setup script for fits2hdf package
"""
from setuptools import setup, find_packages

version = '1.0.0'

# create entry points
# see http://astropy.readthedocs.org/en/latest/development/scripts.html
entry_points = {
    'console_scripts' :
        ['filutil = filterbank.filterbank:cmd_tool',
         'rawutil = filterbank.guppi:cmd_tool',
         'fil2hdf = filterbank.fil2hdf:cmd_tool',
         'gup2hdf = filterbank.gup2hdf:cmd_tool'
     ]
    }

setup(name='filterbank',
      version=version,
      description='Python utilities for Breakthrough Listen SETI observations',
      install_requires=['astropy', 'numpy'],
      url='https://github.com/ucberkeleyseti/filterbank',
      author='Danny Price',
      author_email='dancpr@berkeley.edu',
      license='MIT',
      entry_points=entry_points,
      packages=find_packages(),
      zip_safe=False,
      )
