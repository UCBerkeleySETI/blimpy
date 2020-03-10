[![Build Status](https://travis-ci.org/UCBerkeleySETI/blimpy.svg?branch=master)](https://travis-ci.org/UCBerkeleySETI/blimpy)
[![Documentation Status](https://readthedocs.org/projects/blimpy/badge/?version=latest)](https://blimpy.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/UCBerkeleySETI/blimpy/branch/master/graph/badge.svg)](https://codecov.io/gh/UCBerkeleySETI/blimpy)
 [![JOSS status](http://joss.theoj.org/papers/e58ef21f0a924041bf9438fd75f8aed0/status.svg)](http://joss.theoj.org/papers/e58ef21f0a924041bf9438fd75f8aed0)

## Breakthrough Listen I/O Methods for Python.

### Filterbank + Raw file readers

This repository contains Python 2/3 readers for interacting with [Sigproc filterbank](http://sigproc.sourceforge.net/sigproc.pdf) (.fil), HDF5 (.h5) and [guppi raw](https://baseband.readthedocs.io/en/stable/guppi/) (.raw) files,
as used in the [Breakthrough Listen](https://seti.berkeley.edu) search for intelligent life.


### Installation

The latest release can be installed via pip:

```
pip install blimpy
```

Or, the latest version of the development code can be installed from the github [repo](https://github.com/UCBerkeleySETI/blimpy) and then run `python setup.py install` or `pip install .` (with sudo if required), or by using the following terminal command:

```
pip install https://github.com/UCBerkeleySETI/blimpy/tarball/master
```

To install everything required to run the unit tests, run:

```
pip install -e .[full]
```

You will need `numpy`, `h5py`, `astropy`, `scipy`, and `matplotlib` as dependencies. A `pip install` should pull in numpy, h5py, and astropy, but you may still need to install scipy and matplotlib separately. 
To interact with compressed files, you'll need the `hdf5plugin` package too.

Note that h5py generally needs to be installed in this way:

```
$ pip install --no-binary=h5py h5py
```



### Command line utilities

After installation, some command line utilities will be installed:
* `watutil`, for reading/writing/plotting blimpy filterbank files (either .h5 or .fil format).
* `filutil`, for reading/plotting blimpy filterbank files (.fil format).
* `rawutil`, for plotting data in guppi raw files.
* `fil2h5`, for converting .fil files into .h5 format.
* `h52fil`, for converting .h5 files into .fil format.
* `bldice`, for dicing a smaller frequency region from (either from/to .h5 or .fil).
* `matchfils`, for checking if two .fil files are the same.

Use the `-h` flag to any of the above command line utilities to display their available arguments.

### Reading blimpy filterbank files in .fil or .h5 format

The `blimpy.Waterfall`  provides a Python API for interacting with filterbank data. It supports all BL filterbank data products; see this [example Jupyter notebook](https://github.com/UCBerkeleySETI/blimpy/blob/master/examples/voyager.ipynb) for an overview. 

From the python, ipython or jupiter notebook environments.

```python
from blimpy import Waterfall
fb = Waterfall('/path/to/filterbank.fil')
#fb = Waterfall('/path/to/filterbank.h5') #works the same way
fb.info()
data = fb.data
```

### Reading guppi raw files
The [Guppi Raw format](https://github.com/UCBerkeleySETI/breakthrough/blob/master/doc/RAW-File-Format.md) can be read using the `GuppiRaw` class from `guppi.py`:

```python
from blimpy import GuppiRaw
gr = GuppiRaw('/path/to/guppirawfile.raw')

header, data = gr.read_next_data_block()
```

or

```python
from blimpy import GuppiRaw
gr = GuppiRaw('/path/to/guppirawfile.raw')

for header, data_x, data_y in gr.get_data():
    # process data
```

Note: most users should start analysis with filterbank files, which are smaller in size and have been generated from the guppi raw files.

### Using blimpy inside Docker
The blimpy images are pushed to a public repository after each successful build on Travis.
If you have Docker installed, you can run the following commands to pull our images, which have the environment and dependencies set up for you.

For python3, use:

`docker pull fx196/blimpy:py3_kern_stable`

For python2, use:

`docker pull fx196/blimpy:py2_kern_stable`

Here is a [more complete guide](./docker_guide.md) on using blimpy in Docker.

### Further reading

A detailed overview of the data formats used in Breakthrough Listen can be found in our [data format paper](https://ui.adsabs.harvard.edu/abs/2019arXiv190607391L/abstract). An archive of data files from the Breakthrough Listen program is provided at [seti.berkeley.edu/opendata](http://seti.berkeley.edu/opendata).

### If you have any requests or questions, please lets us know!
