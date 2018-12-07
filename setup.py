"""
setup.py -- setup script for use of packages.
"""
from setuptools import setup, find_packages

__version__ = '1.3.5'

# create entry points
# see http://astropy.readthedocs.org/en/latest/development/scripts.html
entry_points = {
    'console_scripts' : [
        'filutil = blimpy.filterbank:cmd_tool',
        'watutil = blimpy.waterfall:cmd_tool',
        'rawutil = blimpy.guppi:cmd_tool',
        'fil2h5 = blimpy.fil2h5:cmd_tool',
        'h52fil = blimpy.h52fil:cmd_tool',
        'matchfils = blimpy.match_fils:cmd_tool',
        'bldice = blimpy.dice:cmd_tool'
     ]
}

install_requires = [
        'matplotlib<3.0;python_version=="2.7"',
        'matplotlib;python_version>"2.7"',
        'astropy<3.0;python_version=="2.7"',
        'astropy;python_version>"2.7"',
        'numpy',
        'cython',
        'h5py',
        'scipy',
]

extras_require = {
        'full': [
            'bitshuffle',
            'pyslalib',
        ]
}

setup(name='blimpy',
      version=__version__,
      description='Python utilities for Breakthrough Listen SETI observations',
      long_description="Python utilities for Breakthrough Listen SETI observations. It includes data handling, formating, dicing and plotting.",
      platform=['*nix'],
      license='BSD',
      install_requires=install_requires,
      extras_require=extras_require,
      url='https://github.com/ucberkeleyseti/blimpy',
      author='Danny Price, Emilio Enriquez, Griffin Foster, Greg Hellbourg',
      author_email='dancpr@berkeley.edu',
      entry_points=entry_points,
      packages=find_packages(),
      zip_safe=False,
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Natural Language :: English',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 2.7',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Topic :: Scientific/Engineering :: Astronomy',
      ],
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      test_suite="tests",
)
