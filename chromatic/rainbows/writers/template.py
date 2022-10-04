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

# define list of the only things that will show up in imports
__all__ = ["to_abcdefgh"]


def to_abcdefgh(self, filepath):
    """
    Write a Rainbow to a file in the abcdefgh format.

    Parameters
    ----------

    self : Rainbow
        The object to be saved.

    filepath : str
        The path to the file to write.
    """

    # a 1D array of wavelengths (with astropy units of length)
    the_1D_array_of_wavelengths = self.wavelike["wavelength"]

    # a 1D array of times (with astropy units of time)
    the_1D_array_of_times = self.timelike["time"]

    # a 2D (row = wavelength, col = array of fluxes
    the_2D_array_of_fluxes = self.fluxlike["flux"]

    # write out your file, however you like
    write_to_abcdefgh(
        filepath,
        the_1D_array_of_wavelengths,
        the_1D_array_of_times,
        the_2D_array_of_fluxes,
        **some_other_stuff_too_maybe,
    )

    # add some warnings if there's any funny business
    if something_goes_wonky():
        cheerfully_suggest(
            f"""
        Here's a potential problem that the user should know about.
        """
        )
