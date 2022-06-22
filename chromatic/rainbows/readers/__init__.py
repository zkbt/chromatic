from .x1dints import *
from .eureka_txt import *
from .eureka_specdata import *
from .eureka_lcdata import *
from .eureka_channels import *
from .rainbow_npy import *
from .rainbow_FITS import *
from .text import *

# some particular instruments
from .nres import *
from .atoca import *

# some individual-ish folks
from .espinoza import *
from .dossantos import *
from .feinstein import *
from .schlawin import *
from .coloumbe import *


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

    from ...imports import expand_filenames
    from fnmatch import fnmatch
    import glob

    # get all the possible filenames (= expand wildcard)
    filenames = expand_filenames(filepath)
    f = filenames[0]

    # if format='abcdefgh', return the `from_abcdefgh` function
    if format is not None:
        return available_readers[f"from_{format}"]
    # does it look like a .rainbow.npy chromatic file?
    elif fnmatch(f, "*.rainbow.npy"):
        return from_rainbow_npy
    elif fnmatch(f.lower(), "*.rainbow.fits") or fnmatch(f.lower(), "*.rainbow.fit"):
        return from_rainbow_FITS
    # does it look like a .rainbow.npy chromatic file?
    elif fnmatch(f, "*order*.npy"):
        return from_espinoza
    # does it look like a STScI x1dints.fits file?
    elif fnmatch(f, "*x1dints.fits") or fnmatch(f, "*extract_1d.fits"):
        return from_x1dints
    # does it look like an Eureka! S3 text file?
    elif fnmatch(f, "*S3_*_Save.dat") or fnmatch(f, "*S3_*_Save.txt"):
        return from_eureka_SpecData
    # does it look like an Eureka! S3 SpecData hdf5 file?
    elif fnmatch(f, "*S3_*_SpecData.h5"):
        return from_eureka_S3
    # does it look like an Eureka! S4 SpecData hdf5 file?
    elif fnmatch(f, "*S4_*_SpecData.h5"):
        return from_eureka_SpecData
    # does it look like an Eureka! S4 LCData hdf5 file?
    elif fnmatch(f, "*S4_*_LCData.h5"):
        return from_eureka_S4
    elif fnmatch(f, "*S5_*_Table_Save_*.txt"):
        return from_eureka_S5
    elif fnmatch(f, "*.txt") or fnmatch(f, "*.csv"):
        return from_text
    elif fnmatch(f, "*e92-1d.fits") or fnmatch(f, "*e92-1d.fits.fz"):
        return from_nres
    else:
        raise RuntimeError(f"🌈 Failed to guess a good reader for {filenames}.")
