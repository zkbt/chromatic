from .x1dints import *
from .eureka import *
from .rainbow_npy import *
from .rainbow_FITS import *
from .text import *

from .nestor_niriss_numpy import *
from .leodossantos import *
from .adinafeinstein import *
from .everettschlawin import *

# construct a dictionary of available readers
available_readers = {k: globals()[k] for k in globals() if k[0:5] == "from_"}


def guess_reader(filepath, format=None):
    """
    A wrapper to guess the appropriate reader from the filename
    (and possibily an explicitly-set file format string).

    Parameters
    ----------
    filepath : str
        The path to the file, or a glob-friendly string
        pointing to a group of files that should all be
        loaded together (for example, if an exposure was
        split into multiple segments).
    format : str, None
        The file format to use.
    """
    import fnmatch, glob
    from ...imports import expand_filenames

    # get all the possible filenames (= expand wildcard)
    filenames = expand_filenames(filepath)

    # if format='abcdefgh', return the `from_abcdefgh` function
    if format is not None:
        return available_readers[f"from_{format}"]
    # does it look like a .rainbow.npy chromatic file?
    elif fnmatch.fnmatch(filenames[0], "*.rainbow.npy"):
        return from_rainbow_npy
    elif fnmatch.fnmatch(filepath.lower(), "*.rainbow.fits") or fnmatch.fnmatch(
        filepath.lower(), "*.rainbow.fit"
    ):
        return from_rainbow_FITS
    # does it look like a .rainbow.npy chromatic file?
    elif fnmatch.fnmatch(filenames[0], "*order*.npy"):
        return from_nestor_niriss_numpy
    # does it look like a STScI x1dints.fits file?
    elif fnmatch.fnmatch(filenames[0], "*x1dints.fits") or fnmatch.fnmatch(
        filenames[0], "*extract_1d.fits"
    ):
        return from_x1dints
    # does it look like an Eureka! save file?
    elif (
        fnmatch.fnmatch(filenames[0], "*S3_*_Save.dat")
        or fnmatch.fnmatch(filenames[0], "*S3_*_Save.h5")
        or fnmatch.fnmatch(filenames[0], "*S3_*_Save.txt")
    ):
        return from_eureka
    elif fnmatch.fnmatch(filenames[0], "*.txt") or fnmatch.fnmatch(
        filenames[0], "*.csv"
    ):
        return from_text
    else:
        raise RuntimeError(f"🌈 Failed to guess a good reader for {filenames}.")
