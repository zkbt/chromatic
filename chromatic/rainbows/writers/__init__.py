from .rainbow_npy import *
from .rainbow_FITS import *
from .xarray_stellar_spectra import *
from .xarray_raw_light_curves import *
from .xarray_fitted_light_curves import *

from .text import *


# construct a dictionary of available writers
available_writers = {k: globals()[k] for k in globals() if k[0:3] == "to_"}


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
    from fnmatch import fnmatch
    import glob

    # get all the possible filenames (= expand wildcard)
    f = filepath

    # if format='abcdefgh', return the `to_abcdefgh` function
    if format is not None:
        return available_writers[f"to_{format}"]
    # does it look like a .rainbow.npy chromatic file?
    elif fnmatch(f, "*.rainbow.npy"):
        return to_rainbow_npy
    elif fnmatch(f.lower(), "*.rainbow.fits") or fnmatch(f.lower(), "*.rainbow.fit"):
        return to_rainbow_FITS
    # does it look like an ERS-xarray format?
    elif fnmatch(f, "*stellar-spec*.xc"):
        return to_xarray_stellar_spectra
    # does it look like an ERS-xarray format?
    elif fnmatch(f, "*raw-light-curve*.xc"):
        return to_xarray_raw_light_curves
    # does it look like an ERS-xarray format?
    elif fnmatch(f, "*fitted-light-curve*.xc"):
        return to_xarray_fitted_light_curves
    elif fnmatch(f, "*.txt") or fnmatch(filepath, "*.csv"):
        return to_text
    else:
        raise ValueError(
            f"""
        We're having trouble guessing the output format from the filename
        {f}
        Please try specifying a `format=` keyword to your `.save` call.
        """
        )
