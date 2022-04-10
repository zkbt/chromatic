# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_radica"]


def from_radica(rainbow, filepath, order=1):
    """
    Populate a Rainbow from a file in the radica format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load.
    order : int
        The spectral order to be loaded.
    """

    # open the FITS file
    hdu_list = fits.open(filepath)

    # load the header into metadata
    rainbow.metadata = dict(hdu_list["primary"].header)

    # populate a 1D array of wavelengths (with astropy units of length)
    wavelength_2d = hdu_list[f"wo{order}"].data
    rainbow.wavelike["wavelength"] = np.median(wavelength_2d, 0) * u.micron

    # populate a 1D array of times (with astropy units of time)
    # populate a 1D array of times (with astropy units of time)
    times = np.arange(wavelength_2d.shape[0]) * u.minute
    warnings.warn("\nThe times are totally made up!")
    rainbow.timelike["time"] = times

    # populate a 2D (row = wavelength, col = array of fluxes
    rainbow.fluxlike["flux"] = hdu_list[f"fo{order}"].data.T
    rainbow.fluxlike["uncertainty"] = hdu_list[f"eo{order}"].data.T

    message = f"""
    Loading NIRISS order '{order}'. If you want the other order,
    trying `r = Rainbow(..., format='radica', order=2)`
    """
    warnings.warn(message)
