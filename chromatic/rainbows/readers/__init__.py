from .stsci import *
from .eureka import *


def guess_reader(filepath, format=None):
    """
    A wrapper to guess the appropriate from the filename
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

    # get all the possible filenames (= expand wildcard)
    filenames = expand_filenames(filepath)

    # if format='abcdefgh', return the `from_abcdefgh` function
    if format is not None:
        return locals()[f"from_{format}"]
    # does it look like a STScI x1dints.fits file?
    elif fnmatch.fnmatch(filenames[0], "*x1dints.fits"):
        return from_x1dints
    # does it look like an Eureka! save file?
    elif (
        fnmatch.fnmatch(filenames[0], "*S3_*_Save.dat")
        or fnmatch.fnmatch(filenames[0], "*S3_*_Save.h5")
        or fnmatch.fnmatch(filenames[0], "*S3_*_Save.txt")
    ):
        return from_eureka
