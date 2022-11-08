import matplotlib

matplotlib.rcParams["figure.dpi"] = 300
matplotlib.rcParams["figure.constrained_layout.use"] = True
matplotlib.rcParams["figure.facecolor"] = "white"
matplotlib.rcParams["savefig.facecolor"] = "white"
matplotlib.rcParams["savefig.dpi"] = 300

# basics
import numpy as np
import matplotlib.pyplot as plt


import matplotlib.animation as ani
import matplotlib.colors as col
import matplotlib.gridspec as gs

import copy, pkg_resources, os, glob, pickle
from tqdm.auto import tqdm

import warnings, textwrap

from time import time as get_current_seconds


def custom_formatwarning(message, *args, **kwargs):
    return f"🌈🤖 {textwrap.dedent(str(message)).strip().strip()}\n\n"


original_warning_format = warnings.formatwarning


def cheerfully_suggest(*args, **kwargs):
    warnings.formatwarning = custom_formatwarning
    warnings.warn(*args, **kwargs)
    warnings.formatwarning = original_warning_format


"""def custom_formatwarning(message, *args, **kwargs):
    # ignore everything except the message
    return f"\n🌈 Warning: {textwrap.dedent(str(message))}"

warnings.formatwarning = custom_formatwarning
"""
# astropy
from astropy.io import ascii, fits
from astropy.table import Table, QTable
from astropy.time import Time
from astropy.stats import sigma_clip, median_absolute_deviation, mad_std

# import astropy.units as u
import astropy.constants as con
from astropy.visualization import quantity_support, simple_norm

# slightly fancier visualization tools
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from scipy.interpolate import interp1d

# For converting Rainbows to pandas dataframe
import pandas as pd

from .units import *

# For remove astrophysical signal
from scipy.signal import savgol_filter, butter, filtfilt
from scipy.signal import medfilt, convolve2d
from scipy.ndimage import median_filter

# define a directory where we can put any necessary data files
data_directory = pkg_resources.resource_filename("chromatic", "data")


def is_being_run_from_jupyter():
    try:
        return get_ipython().__class__.__name__ == "ZMQInteractiveShell"
    except NameError:
        return False


from astropy.utils.data import download_file


def download_file_with_warning(*args, **kwargs):
    if is_being_run_from_jupyter():
        cheerfully_suggest(
            """
        The progress bar for this download is not being shown
        because `astropy.utils.data.download_file` is being run
        from a jupyter notebook instead of from the terminal.
        """
        )
    kwargs["show_progress"] = True
    return download_file(*args, **kwargs)


def expand_filenames(filepath):
    """
    A wrapper to expand a string or list into a list of filenames.
    """
    if type(filepath) == list:
        filenames = filepath
    elif "*" in filepath:
        filenames = np.sort(glob.glob(filepath))
    else:
        filenames = [filepath]
    return sorted(filenames)


def name2color(name):
    """
    Return the 3-element RGB array of a given color name.

    Parameters
    ----------
    name : str
        The name of a color

    Returns
    -------
    rgb : tuple
        3-element RGB color, with numbers from 0.0 to 1.0
    """

    # give a friendly warning if the color name can't be found
    try:
        color_hex = col.cnames[name]
        return col.hex2color(color_hex)
    except KeyError:
        cheerfully_suggest(f"The color {name} can't be found. (Returning black.)")
        return (0.0, 0.0, 0.0)


def one2another(bottom="white", top="red", alpha_bottom=1.0, alpha_top=1.0, N=256):
    """
    Create a cmap that goes smoothly (linearly in RGBA) from "bottom" to "top".

    Parameters
    ----------
    bottom : str
        Name of a color for the bottom of cmap (0.0)
    top : str
        Name of a color for the top of the cmap (1.0)
    alpha_bottom : float
        Opacity at the bottom of the cmap
    alpha_top : float
        Opacitiy at the top of the cmap
    N : int
        The number of levels in the listed color map

    Returns
    -------
    cmap : Colormap
        A color map that goes linearly from the
        bottom to top color (and alpha).
    """

    # get the RGB values of the bottom and top of the cmap
    rgb_bottom, rgb_top = name2color(bottom), name2color(top)

    # create linear gradients for all four RGBA channels
    r = np.linspace(rgb_bottom[0], rgb_top[0], N)
    g = np.linspace(rgb_bottom[1], rgb_top[1], N)
    b = np.linspace(rgb_bottom[2], rgb_top[2], N)
    a = np.linspace(alpha_bottom, alpha_top, N)

    # create (N,4) array + populate a listed colormap
    colors = np.transpose(np.vstack([r, g, b, a]))
    cmap = col.ListedColormap(colors, name="{bottom}2{top}".format(**locals()))

    # return the colormap
    return cmap


def remove_unit(x):
    """
    Quick wrapper to remove the unit from a quantity,
    but not complain if it doesn't have one.
    """
    try:
        return x.value
    except AttributeError:
        return x
