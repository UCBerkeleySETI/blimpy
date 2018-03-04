"""
setup.py -- setup script for use of packages.
"""
from setuptools import setup, find_packages

__version__ = '1.1.7'

# create entry points
# see http://astropy.readthedocs.org/en/latest/development/scripts.html
entry_points = {
    'console_scripts' :
        ['filutil = blimpy.filterbank:cmd_tool',
         'watutil = blimpy.waterfall:cmd_tool',
         'rawutil = blimpy.guppi:cmd_tool',
         'fil2h5 = blimpy.fil2h5:cmd_tool',
         'h52fil = blimpy.h52fil:cmd_tool',
         'matchfils = blimpy.match_fils:cmd_tool'
     ]
    }

setup(name='blimpy',
    version = __version__,
    description = 'Python utilities for Breakthrough Listen SETI observations',
    long_description = "Python utilities for Breakthrough Listen SETI observations",
    platforms = ['*nix'],
    license = 'MIT',
    install_requires = ['astropy', 'numpy', 'cython', 'h5py'],
    url='https://github.com/ucberkeleyseti/blimpy',
    author='Danny Price, Emilio Enriquez, Griffin Foster',
    author_email='dancpr@berkeley.edu',
    entry_points=entry_points,
    packages=find_packages(),
    zip_safe=False,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)

