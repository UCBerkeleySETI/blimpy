"""
setup.py -- setup script for use of packages.
"""
from setuptools import setup, find_packages

__version__ = '2.0.0'

with open("README.md", "r") as fh:
    long_description = fh.read()

# create entry points
# see http://astropy.readthedocs.org/en/latest/development/scripts.html
entry_points = {
    'console_scripts' : [
        'filutil = blimpy.filterbank:cmd_tool',
        'watutil = blimpy.waterfall:cmd_tool',
        'rawutil = blimpy.guppi:cmd_tool',
        'fil2h5 = blimpy.fil2h5:cmd_tool',
        'h52fil = blimpy.h52fil:cmd_tool',
        'bl_scrunch = blimpy.bl_scrunch:cmd_tool',
        'matchfils = blimpy.match_fils:cmd_tool',
        'bldice = blimpy.dice:cmd_tool'
     ]
}

install_requires = [
        'setuptools<40.0;python_version=="2.7"',
        'matplotlib<3.0;python_version=="2.7"',
        'matplotlib;python_version>"2.7"',
        'astropy<3.0;python_version=="2.7"',
        'astropy;python_version>"2.7"',
        'numpy',
        'cython',
        'h5py',
        'scipy',
        'six',
        'hdf5plugin',
        'pandas'
]

extras_require = {
        'full': [
            'pyslalib',
        ]
}

setup(name='blimpy',
      version=__version__,
      description='Python utilities for Breakthrough Listen SETI observations',
      long_description=long_description,
      long_description_content_type='text/markdown',
      license='BSD',
      install_requires=install_requires,
      extras_require=extras_require,
      url='https://github.com/ucberkeleyseti/blimpy',
      author='Danny Price, Emilio Enriquez, Yuhong Chen, Mark Siebert, and BL contributors',
      author_email='dancpr@berkeley.edu',
      entry_points=entry_points,
      packages=find_packages(),
      include_package_data=True,
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
      test_suite="blimpytests",
)
