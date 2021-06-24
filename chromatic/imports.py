# basics
import numpy as np
import matplotlib.pyplot as plt
import copy, pkg_resources, os, glob, warnings
from tqdm import tqdm
import copy

# astropy
from astropy.io import ascii, fits
import astropy.units as u
import astropy.constants as con
from astropy.visualization import quantity_support

# slightly fancier visualization tools
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from scipy.interpolate import interp1d


# define a driectory where we can put any necessary data files
data_directory = pkg_resources.resource_filename("chromatic", "data")
