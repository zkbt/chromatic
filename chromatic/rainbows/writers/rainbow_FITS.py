"""
Define a writer for chromatic .rainbow.FITS files.
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["to_rainbow_FITS"]


def to_rainbow_FITS(rainbow, filepath, overwrite=None):
    """
    Write a Rainbow to a FITS file.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be saved.

    filepath : str
        The path to the file to write.
    """

    # create a header for the metadata
    header = fits.Header()
    for k in rainbow.metadata:
        header[k] = rainbow.metadata[k]
    primary_hdu = fits.PrimaryHDU(header=header)

    # create extensions for the three other core dictionaries
    flux_hdu = fits.BinTableHDU(Table(rainbow.fluxlike), name="fluxlike")
    wave_hdu = fits.BinTableHDU(Table(rainbow.wavelike), name="wavelike")
    time_hdu = fits.BinTableHDU(Table(rainbow.timelike), name="timelike")

    hdu_list = fits.HDUList([primary_hdu, flux_hdu, wave_hdu, time_hdu])
    hdu_list.writeto(filepath, overwrite=overwrite)
