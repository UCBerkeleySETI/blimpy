## Filterbank + RAW file readers

This repository contains basic python readers for interacting with filterbank (.fil) and guppi raw (.raw) files generated
by the Breakthrough Listen instruments.

### Reading filterbank files
The [Sigproc user guide](http://sigproc.sourceforge.net/sigproc.pdf) gives details of the filterbank format. The `filterbank.py` script provides a Python API for interacting with filterbank data; see this [example Jupyter notebook](https://github.com/UCBerkeleySETI/breakthrough/blob/master/GBT/voyager/voyager.ipynb) for an overview.

### Reading guppi raw files
The [Guppi Raw format](https://github.com/UCBerkeleySETI/breakthrough/blob/master/doc/RAW-File-Format.md) can be read using the `GuppiRaw` class from `guppi.py`:

```python
from guppi import GuppiRaw
r = GuppiRaw('/path/to/guppirawfile.raw')

header, data = r.read_next_data_block()
```

Most users should start analysis with filterbank files, which are smaller in size and have been generated from the guppi raw files.
