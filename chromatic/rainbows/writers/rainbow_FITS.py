"""
Define a writer for chromatic .rainbow.FITS files.
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["to_rainbow_FITS"]


def to_rainbow_FITS(self, filepath, overwrite=True):
    """
    Write a Rainbow to a FITS file.

    Parameters
    ----------

    self : Rainbow
        The object to be saved.

    filepath : str
        The path to the file to write.
    """

    # create a header for the metadata
    header = fits.Header()
    header["comment"] = "Data are stored in three main extensions"
    header["comment"] = " [1]=['FLUXLIKE'] contains (nwave, ntime)-shaped arrays"
    header["comment"] = " [2]=['WAVELIKE'] contains (nwave)-shaped arrays"
    header["comment"] = " [3]=['TIMELIKE'] contains (ntime)-shaped arrays"
    header["comment"] = "The primary extension header contains some metadata."

    for k in self.metadata:
        try:
            header[k] = self.metadata[k]
        except ValueError:
            cheerfully_suggest(f"metadata item '{k}' cannot be saved to FITS header")

    primary_hdu = fits.PrimaryHDU(header=header)

    # create extensions for the three other core dictionaries
    flux_hdu = fits.BinTableHDU(Table(self.fluxlike), name="fluxlike")
    wave_hdu = fits.BinTableHDU(Table(self.wavelike), name="wavelike")
    time_hdu = fits.BinTableHDU(Table(self.timelike), name="timelike")

    hdu_list = fits.HDUList([primary_hdu, flux_hdu, wave_hdu, time_hdu])
    hdu_list.writeto(filepath, overwrite=overwrite)
