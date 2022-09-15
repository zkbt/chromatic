# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_eureka_channels", "from_eureka_S5"]


def from_eureka_channels(rainbow, filepath):
    """
    Populate a Rainbow from Eureka! S5 modeling outputs,
    which are wavelength-binned light curves with models.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.

    filepath : str
        The path to the files to load.

    """
    # get list of files
    filenames = expand_filenames(filepath)

    # load them into tables
    tables = [ascii.read(f) for f in filenames]

    # load the time and wavelength arrays
    rainbow.timelike["time"] = tables[0]["time"] * u.day
    rainbow.wavelike["wavelength"] = [t["wavelength"][0] for t in tables] * u.micron
    rainbow.wavelike["wavelength_lower"] = [
        (t["wavelength"] - t["bin_width"])[0] for t in tables
    ] * u.micron
    rainbow.wavelike["wavelength_upper"] = [
        (t["wavelength"] + t["bin_width"])[0] for t in tables
    ] * u.micron
    keys_used = ["time", "wavelength", "bin_width"]

    # populate the fluxlike arrays
    for k in ["lcdata", "lcerr", "model", "residuals"]:
        rainbow.fluxlike[k] = np.array([t[k] for t in tables])
        keys_used.append(k)

    for k in tables[0].colnames:
        if k != keys_used:
            rainbow.fluxlike[k] = np.array([t[k] for t in tables])

    rainbow.fluxlike["flux"] = rainbow.lcdata * 1
    rainbow.fluxlike["uncertainty"] = rainbow.lcerr * 1


from_eureka_S5 = from_eureka_channels
