"""
Define a writer for chromatic .rainbow.npy files.
"""


# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["to_text"]


def to_text(rainbow, filepath, overwrite=None, group_by="wavelength"):
    """
    Write a Rainbow to a file in the text format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be saved.

    filepath : str
        The path to the file to write.
    """

    # a 1D array of wavelengths (with astropy units of length)
    the_1D_array_of_wavelengths = rainbow.wavelike["wavelength"]

    # a 1D array of times (with astropy units of time)
    the_1D_array_of_times = rainbow.timelike["time"]

    # a 2D (row = wavelength, col = array of fluxes
    the_2D_array_of_fluxes = rainbow.fluxlike["flux"]

    # write out your file, however you like
    w, t = np.meshgrid(
        rainbow.wavelength.to("micron"), rainbow.time.to("day"), indexing="ij"
    )

    def make_into_columns(x):
        if group_by == "wavelength":
            return x.flatten()
        if group_by == "time":
            return x.T.flatten()

    # create a table
    table = Table(
        dict(wavelength=make_into_columns(w), time=make_into_columns(t)),
        meta=rainbow.metadata,
    )
    table["flux"] = make_into_columns(rainbow.flux)
    try:
        table["uncertainty"] = make_into_columns(rainbow.uncertainty)
    except KeyError:
        pass

    other_keys = list(rainbow.fluxlike.keys())
    other_keys.remove("flux")
    other_keys.remove("uncertainty")

    for k in other_keys:
        table[k] = make_into_columns(rainbow.fluxlike[k])

    table.write(filepath, format="ascii.ecsv", overwrite=overwrite)
