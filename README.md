[![Build Status](https://travis-ci.org/UCBerkeleySETI/blimpy.svg?branch=master)](https://travis-ci.org/UCBerkeleySETI/blimpy)
[![Coverage Status](https://coveralls.io/repos/github/telegraphic/blimpy/badge.svg?branch=master)](https://coveralls.io/github/telegraphic/blimpy?branch=master)

## Breakthrough Listen I/O Methods for Python.

### Filterbank + Raw file readers

This repository contains Python 2 readers for interacting with Sigproc filterbank (.fil), HDF5 (.h5) and guppi raw (.raw) files generated
by the Breakthrough Listen instruments.

### Installation

The latest stable release (1.2.0) can be installed via pip:

```
pip install blimpy
```

Or, the latest version of the development code can be installed from the github [repo](https://github.com/UCBerkeleySETI/blimpy) and then run `python setup.py install` (with sudo if required) or by using the following terminal command:

```
pip install https://github.com/UCBerkeleySETI/blimpy/tarball/master
```

You will need numpy, h5py astropy, scipy, and matplotlib as dependencies. A `pip install` should pull in numpy, h5py, and astropy, but you may still need to install scipy and matplotlib separately.

```
Version 1.2.0
Now more support/compatibility for Python 2 and 3.
Version 1.1.9
Now supporting dicing (bldice).
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

The `blimpy.Waterfall`  provides a Python API for interacting with filterbank data. It supports all BL filterbank data products; see this [example Jupyter notebook](https://github.com/UCBerkeleySETI/breakthrough/blob/master/GBT/voyager/voyager.ipynb) for an overview.
The [Sigproc user guide](http://sigproc.sourceforge.net/sigproc.pdf) gives details of the filterbank (.fil) format. This method is the most up-to-date and supercedes the depreciated `blimpy.Filterbank` which can only import .fil files.

From the python, ipython or jupiter notebook environments.
```
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
r = GuppiRaw('/path/to/guppirawfile.raw')

header, data = r.read_next_data_block()
```

Note: most users should start analysis with filterbank files, which are smaller in size and have been generated from the guppi raw files.


### This readme is far from complete. If you have any request/questions please lets us know!
