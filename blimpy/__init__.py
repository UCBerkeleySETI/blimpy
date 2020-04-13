from __future__ import absolute_import


try:
    from . import waterfall
    from .waterfall import Waterfall
    from .guppi import GuppiRaw
    from . import utils
except ImportError:
    print("Warning: Cannot import main utilities")

try:
    from . import fil2h5
    from . import h52fil
    from . import bl_scrunch
except:
    print("Warning: Cannot import HDF5 utilities")

try:
    from blimpy.io import file_wrapper
    from . import match_fils
except:
    pass

from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution('blimpy').version
except DistributionNotFound:
    __version__ = '0.0.0 - please install via pip/setup.py'
