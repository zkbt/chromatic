from .rainbow_npy import *
from .rainbow_FITS import *
from .text import *


# construct a dictionary of available writers
available_writers = {k: globals()[k] for k in globals() if k[0:5] == "to_"}


def guess_writer(filepath, format=None):
    """
    A wrapper to guess the appropriate writer from the filename
    (and possibily an explicitly-set file format string).

    Parameters
    ----------
    filepath : str
        The path to the file to be written.
    format : str, None
        The file format to use.
    """
    import fnmatch, glob

    # if format='abcdefgh', return the `to_abcdefgh` function
    if format is not None:
        return available_writers[f"to_{format}"]
    # does it look like a .rainbow.npy chromatic file?
    elif fnmatch.fnmatch(filepath, "*.rainbow.npy"):
        return to_rainbow_npy
    elif fnmatch.fnmatch(filepath.lower(), "*.rainbow.fits") or fnmatch.fnmatch(
        filepath.lower(), "*.rainbow.fit"
    ):
        return to_rainbow_FITS
    elif fnmatch.fnmatch(filepath, "*.txt") or fnmatch.fnmatch(filepath, "*.csv"):
        return to_text
