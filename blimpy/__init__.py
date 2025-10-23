from __future__ import absolute_import

try:
    from . import waterfall
    from .waterfall import Waterfall
    from .guppi import GuppiRaw
    from . import utils
    from . import fil2h5
    from . import h52fil
    from . import h5diag
    from . import bl_scrunch
    from . import calcload
    from . import rawhdr
    from . import stax
    from . import stix
    from . import dsamp
    from blimpy.io import file_wrapper
except:
    print("Warning: At least one utility could not be imported!")

try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # Fallback for Python < 3.8
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version('blimpy')
except PackageNotFoundError:
    __version__ = '0.0.0 - please install via pip/setup.py'
