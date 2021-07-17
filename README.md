[![Build Status](https://github.com/UCBerkeleySETI/blimpy/workflows/Test%20Blimpy/badge.svg)](https://github.com/UCBerkeleySETI/blimpy/actions)
[![Documentation Status](https://readthedocs.org/projects/blimpy/badge/?version=latest)](https://blimpy.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/UCBerkeleySETI/blimpy/branch/master/graph/badge.svg)](https://codecov.io/gh/UCBerkeleySETI/blimpy)
 [![JOSS status](http://joss.theoj.org/papers/e58ef21f0a924041bf9438fd75f8aed0/status.svg)](http://joss.theoj.org/papers/e58ef21f0a924041bf9438fd75f8aed0)

## Breakthrough Listen I/O Methods for Python.

### Filterbank + Raw file readers

This repository contains Python 2/3 readers for interacting with [Sigproc filterbank](http://sigproc.sourceforge.net/sigproc.pdf) (.fil), HDF5 (.h5) and [guppi raw](https://baseband.readthedocs.io/en/stable/guppi/) (.raw) files,
as used in the [Breakthrough Listen](https://seti.berkeley.edu) search for intelligent life.


### Installation

#### System Dependencies
Sometimes the `pip` installation can fail if a system dependency is not installed. To fix this, make sure you have `curl` and install the required system dependencies with the command bellow:

##### Debian/Ubuntu
```
curl https://raw.githubusercontent.com/UCBerkeleySETI/blimpy/master/dependencies.txt | xargs -n 1 sudo apt install --no-install-recommends -y
```

#### Manual Installation

The latest release can be installed via pip directly from this repository:

```
python3 -m pip install -U git+https://github.com/UCBerkeleySETI/blimpy
```

Or, the latest version of the development code can be installed from the github [repo](https://github.com/UCBerkeleySETI/blimpy) and then run `python setup.py install` or `pip install .` (with sudo if required), or by using the following terminal command:

```
python3 -m pip install -U https://github.com/UCBerkeleySETI/blimpy/tarball/master
```

To install everything required to run the unit tests, run:

```
python3 -m pip install -e .[full]
```

You will need `numpy`, `h5py`, `astropy`, `scipy`, and `matplotlib` as dependencies. A `pip install` should pull in numpy, h5py, and astropy, but you may still need to install scipy and matplotlib separately.
To interact with compressed files, you'll need the `hdf5plugin` package too.

Note that h5py generally needs to be installed in this way:

```
$ python3 -m pip install --no-binary=h5py h5py
```

### Command line utilities

After installation, some command line utilities will be installed:
* `watutil`, Read/write/plot an .h5 file or a .fil file.
* `rawutil`, Plot data in a guppi raw file.
* `fil2h5`, Convert a .fil file into .h5 format.
* `h52fil`, Convert an .h5 file into .fil format.
* `bldice`, Dice a smaller frequency region from (either from/to .h5 or .fil).
* `matchfils`, Check if two .fil files are the same.
* `calcload`, Calculate the Waterfall max_load value needed to load the data array for a given file.
* `rawhdr`, Display the header fields of a raw guppi file.
* `stax`, For a collection of .h5 or .fil files, create a stack of waterfall plots as a PNG file.

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

`docker pull fx196/blimpy:py3_kern_stable`

Here is a [more complete guide](./docker_guide.md) on using blimpy in Docker.

### Further reading

A detailed overview of the data formats used in Breakthrough Listen can be found in our [data format paper](https://ui.adsabs.harvard.edu/abs/2019arXiv190607391L/abstract). An archive of data files from the Breakthrough Listen program is provided at [seti.berkeley.edu/opendata](http://seti.berkeley.edu/opendata).

### If you have any requests or questions, please lets us know!
