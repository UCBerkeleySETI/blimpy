try:
    from pyslalib import slalib as s
    HAS_SLALIB = True
except ImportError:
    HAS_SLALIB = False

import numpy as np
from astropy.coordinates import Angle

# Telescope coordinates (needed for LSR calc)
parkes_coords = (-32.998370, 148.263659,  324.00)
gbt_coords    = (38.4331294, 79.8398397, 824.36)
