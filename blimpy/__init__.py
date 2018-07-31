from __future__ import absolute_import

from .filterbank import Filterbank, read_header, fix_header
from .guppi import GuppiRaw
from . import utils

try:
    from . import fil2h5
    from . import h52fil
    from . import gup2hdf
except:
    print("Warning: Cannot import HDF5 utilities")

from . import file_wrapper
from . import match_fils

try:
    from . import waterfall
    from .waterfall import Waterfall
except ImportError:
    pass

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution('blimpy').version
except DistributionNotFound:
    __version__ = '0.0.0 - please install via pip/setup.py'