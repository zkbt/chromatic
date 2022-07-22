# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_feinstein"]


def from_feinstein(rainbow, filepath):
    """
    Populate a Rainbow from a file in the feinstein format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load.
    """

    wavelength, spectra, err, time = np.load(filepath, allow_pickle=True)

    rainbow.wavelike["wavelength"] = wavelength * u.micron * 1

    # populate a 1D array of times (with astropy units of time)
    times = time * u.day
    rainbow.timelike["time"] = times * 1

    # populate a 2D (row = wavelength, col = array of fluxes
    flux = np.zeros((len(wavelength), len(times)))
    uncertainty = np.zeros_like(flux)

    rainbow.fluxlike["flux"] = spectra.T * 1
    rainbow.fluxlike["uncertainty"] = err.T * 1
