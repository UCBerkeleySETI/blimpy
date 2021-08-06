from __future__ import absolute_import

try:
    from . import waterfall
    from .waterfall import Waterfall
    from .guppi import GuppiRaw
    from . import utils
    from . import fil2h5
    from . import h52fil
    from . import bl_scrunch
    from . import calcload
    from . import rawhdr
    from . import stax
    from . import stix
    from . import match_fils
    from blimpy.io import file_wrapper
except:
    print("Warning: At least one utility could not be imported!")

from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution('blimpy').version
except DistributionNotFound:
    __version__ = '0.0.0 - please install via pip/setup.py'
