# basics
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.colors as col
import matplotlib.gridspec as gs

import copy, pkg_resources, os, glob, fnmatch, pickle
from tqdm import tqdm

import warnings, textwrap


def custom_formatwarning(message, *args, **kwargs):
    # ignore everything except the message
    return f"\n🌈 Warning: {textwrap.dedent(str(message))}"


warnings.formatwarning = custom_formatwarning

# astropy
from astropy.io import ascii, fits
import astropy.units as u
import astropy.constants as con
from astropy.visualization import quantity_support

# slightly fancier visualization tools
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from scipy.interpolate import interp1d

# For modelling transits.
import batman

# define a driectory where we can put any necessary data files
data_directory = pkg_resources.resource_filename("chromatic", "data")


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
    return filenames
