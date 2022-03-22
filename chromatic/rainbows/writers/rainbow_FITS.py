"""
This module serves as a template for creating a new Rainbow
writer. If you want to add the ability to writer chromatic
light curves to a new kind of file format, a good process
would be to do something like the following

    1. Copy this `template.py` file into a new file in the
    `writers/` directory, ideally with a name that's easy
    to recognize, such as `writers/abcdefgh.py` (assuming
    `abcdefgh` is the name of your format)

    2. Start by finding and replacing `abcdefgh` in this
    template with the name of your format.

    3. Edit the `to_abcdefgh` function so that it will
    write out a Rainbow to whatever format you want.

    4. Edit the `writers/__init__.py` file to import your
    `to_abcdefgh` function to be accessible when people
    try to write Rainbows out to files. Add an `elif` statement
    to the `guess_writer` function that will help guess which
    writer to use from some aspect(s) of the filename.

    (This `guess_writer` function also accepts a `format=`
    keyword that allows the user to explicitly specify that
    the abcdefgh writer should be used.)

    5. Submit a pull request to the github repository for
    this package, so that other folks can use your handy
    new writer too!
"""

# import the general list of packages
from ...imports import *

"""
Define a writer for chromatic .rainbow.FITS files.
"""


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
