# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_everett"]


def from_everett(rainbow, filepath):
    """
    Populate a Rainbow from a file in Everett Schlawin's tshirt format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load.
    """

    # read in your file, however you like
    hdu_list = fits.open(filepath)

    # populate a 1D array of wavelengths (with astropy units of length)
    rainbow.wavelike["indices"] = hdu_list["disp indices"].data
    rainbow.wavelike["wavelength"] = hdu_list["wavelength"].data * u.micron

    # populate a 1D array of times (with astropy units of time)
    rainbow.timelike["time"] = hdu_list["time"].data * u.day

    # populate a 2D (row = wavelength, col = array of fluxes
    for k in [
        "optimal spec",
        "opt spec err",
        "sum spec",
        "sum spec err",
        "background spec",
        "refpix",
    ]:
        rainbow.fluxlike[k] = hdu_list[k].data.T.squeeze()
    rainbow.fluxlike["flux"] = rainbow.fluxlike["optimal spec"]
    rainbow.fluxlike["uncertainty"] = rainbow.fluxlike["opt spec err"]

    # kludgily pull all extension headers into the metadata
    for k in hdu_list:
        rainbow.metadata[f"{k}-header"] = hdu_list[k].header
