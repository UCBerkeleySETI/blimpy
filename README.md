[![Build Status](https://travis-ci.org/UCBerkeleySETI/blimpy.svg?branch=master)](https://travis-ci.org/UCBerkeleySETI/blimpy)

## Filterbank + RAW file readers

This repository contains basic python readers for interacting with filterbank (.fil) and guppi raw (.raw) files generated
by the Breakthrough Listen instruments.

### Installation

To install, download this repository and then run `python setup.py install` (with sudo if required). You will need numpy and astropy as dependencies.

### Command line utilities

After installation, two command line utilities will be installed:
* `filutil`, for plotting and reading filterbank information.
* `filutil2`, for writing h5 files (but plotting also available).
* `rawutil`, for plotting data in guppi raw files.

Use the `-h` flag to display command line arguments.

### Reading filterbank files
The [Sigproc user guide](http://sigproc.sourceforge.net/sigproc.pdf) gives details of the filterbank format. The `filterbank.py` script provides a Python API for interacting with filterbank data; see this [example Jupyter notebook](https://github.com/UCBerkeleySETI/breakthrough/blob/master/GBT/voyager/voyager.ipynb) for an overview.

```python
from blimpy import Filterbank
fb = Filterbank('/path/to/filterbank.fil')
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
