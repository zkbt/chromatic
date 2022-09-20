"""
Define a reader for chromatic .rainbow.FITS files.
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_rainbow_FITS"]


def from_rainbow_FITS(rainbow, filepath):
    """
    Populate a Rainbow from a file in the rainbow_FITS format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.

    filepath : str
        The path to the file to load.
    """

    # open the FITS file
    hdu_list = fits.open(filepath)

    # load the header into metadata
    h = hdu_list["primary"].header
    for k in h:
        if k.lower() in [
            "name",
            "wscale",
            "tscale",
            "signal_to_noise",
            "history",
        ]:  # KLUDGE!
            rainbow.metadata[k.lower()] = h[k]
        elif k not in ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]:
            rainbow.metadata[k] = h[k]

    # load into the three core dictionaries
    for e in ["fluxlike", "wavelike", "timelike"]:
        table = Table.read(hdu_list[e])
        for k in table.colnames:
            vars(rainbow)[e][k] = table[k].quantity * 1
