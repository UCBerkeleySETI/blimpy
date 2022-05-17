"""
Blimpy Plotting Configuration
This file is imported by the other plotting source files and blimpy/waterfall.py.
matplotlib backends info:
https://matplotlib.org/3.5.0/users/explain/backends.html#:~:text=By%20default%2C%20Matplotlib%20should%20automatically,to%20worry%20about%20the%20backend.
"""
import os
import numpy as np
import matplotlib

# Define plt for caller.
import matplotlib.pyplot as plt
plt.rcParams['axes.formatter.useoffset'] = False

# Define NullFormatter for caller.
from matplotlib.ticker import NullFormatter

#Define some constants for caller.
MAX_PLT_POINTS      = 65536                  # Max number of points in matplotlib plot
MAX_IMSHOW_POINTS   = (8192, 4096)           # Max number of points in imshow plot


def ok_to_show():
    """
    Tell caller if the DISPLAY environment variable is set
    and therefore if plt.show() can be executed.
    Parameters
    ----------
    None.
    Returns
    -------
    bool
        Can plt.show() be executed (True/False)?
    """
    display = os.environ.get("DISPLAY", "empty")
    if display == "empty":
        print("blimpy plotting config.py setup_plotting_backend: DISPLAY is *empty*")
        return False
    print(f"blimpy plotting config.py setup_plotting_backend: DISPLAY is {display}")
    return True


def print_plotting_backend(arg_context):
    """ Show which matplotlib backend is in use."""
    ok_to_show()
    print(f"blimpy plotting config.py ({arg_context}): matplotlib backend is {matplotlib.get_backend()}")


def get_mpl_backend():
    return matplotlib.get_backend()


def set_mpl_backend(backend):
    matplotlib.use(backend)

