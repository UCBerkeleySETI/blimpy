"""
setup.py -- setup script for use of packages.
"""
from setuptools import setup, find_packages

__version__ = '2.1.1'

with open("README.md", "r") as fh:
    long_description = fh.read()

# create entry points
# see http://astropy.readthedocs.org/en/latest/development/scripts.html
entry_points = {
    'console_scripts' : [
        'bldice = blimpy.dice:cmd_tool',
        'bl_scrunch = blimpy.bl_scrunch:cmd_tool',
        'calcload = blimpy.calcload:cmd_tool',
        'dsamp = blimpy.dsamp:cmd_tool',
        'fil2h5 = blimpy.fil2h5:cmd_tool',
        'h52fil = blimpy.h52fil:cmd_tool',
        'h5diag = blimpy.h5diag:cmd_tool',
        'matchfils = blimpy.match_fils:cmd_tool',
        'peek = blimpy.peek:cmd_tool',
        'rawhdr = blimpy.rawhdr:cmd_tool',
        'rawutil = blimpy.guppi:cmd_tool',
        'srcname = blimpy.srcname:cmd_tool',
        'stax = blimpy.stax:cmd_tool',
        'stix = blimpy.stix:cmd_tool',
        'watutil = blimpy.waterfall:cmd_tool',
     ]
}

with open("requirements.txt", "r") as fh:
    install_requires = fh.readlines()

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
          'Programming Language :: Python :: 3.7',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Topic :: Scientific/Engineering :: Astronomy',
      ],
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'pyslalib'],
      test_suite="blimpytests",
)
