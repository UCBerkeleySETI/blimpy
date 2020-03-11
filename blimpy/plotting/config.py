import os
import numpy as np

# Check if $DISPLAY is set (for handling plotting on remote machines with no X-forwarding)
import matplotlib
HAS_DISPLAY = True
try:
    from tkinter import TclError
except:
    HAS_DISPLAY = False

if 'DISPLAY' in os.environ.keys() and HAS_DISPLAY:
    try:
        import pylab as plt
    except TclError:
        HAS_DISPLAY = False

if not HAS_DISPLAY:
    matplotlib.use('Agg')
    import pylab as plt

from matplotlib.ticker import NullFormatter

plt.rcParams['axes.formatter.useoffset'] = False


MAX_PLT_POINTS      = 65536                  # Max number of points in matplotlib plot
MAX_IMSHOW_POINTS   = (8192, 4096)           # Max number of points in imshow plot