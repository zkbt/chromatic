"""
Define a reader for Néstor's NIRISS extracted spectra .npy files.

This was kludged together to read the Stage 2 outputs from
https://stsci.app.box.com/s/tyg3qqd85601gkbw5koowrx0obekeg0m/folder/154382588636
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_espinoza"]


def from_espinoza(rainbow, filepath):
    """
    Populate a Rainbow from a file in the espinoza format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load. It should be a glob-friendly
        string like `.../*_order1.npy` that points to both a
        `spectra_order1.npy` and a `wavelengths_order1.npy` file.
    """

    # read in two files, one for the flux, one for the wavelength
    spectra_filename, wavelengths_filename = sorted(glob.glob(filepath))
    assert ("spectra") in spectra_filename
    assert ("wavelength") in wavelengths_filename
    spectra = np.load(spectra_filename)
    wavelengths = np.load(wavelengths_filename)

    # populate a 1D array of wavelengths (with astropy units of length)
    rainbow.wavelike["wavelength"] = wavelengths * u.micron * 1

    # populate a 1D array of times (with astropy units of time)
    times = np.arange(spectra.shape[0]) * u.minute
    cheerfully_suggest("The times are totally made up!")
    rainbow.timelike["time"] = times * 1

    # populate a 2D (row = wavelength, col = array of fluxes
    rainbow.fluxlike["flux"] = spectra[:, 0, :].T * 1
    rainbow.fluxlike["uncertainty"] = spectra[:, 1, :].T * 1
