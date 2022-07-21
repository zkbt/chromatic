from ...imports import *

__all__ = ["from_eureka_S3_txt"]


def eureadka_txt(filename):
    """
    Read eureka's concatenated results table
    and parse it into a Rainbow-friendly format

    Parameters
    ----------
    filename : str
        The path to the file.
    """

    # load the data
    data = ascii.read(filename)

    # figure out a time array
    for time_key in ["time", "bjdtdb"]:
        try:
            t = np.unique(data[time_key])
            break
        except KeyError:
            pass
    timelike = {}
    timelike["time"] = t * u.day * 1

    # figure out a wavelength array
    w = np.unique(data["wave_1d"])
    wavelike = {}
    wavelike["wavelength"] = w * u.micron * 1

    # populate the fluxlike quantities
    fluxlike = {}
    i_time = np.arange(len(t))
    # loop through wavelengths, populating all times for each
    for i_wavelength in tqdm(range(len(w)), leave=False):

        for k in data.colnames[2:]:
            # if an array for this key doesn't exist, create it
            if k not in fluxlike:
                fluxlike[k] = np.zeros((len(w), len(t)))

            # figure out all indices for this wavelengh
            indices_for_this_wavelength = i_wavelength + i_time * len(w)
            fluxlike[k][i_wavelength, i_time] = data[k][indices_for_this_wavelength]

    fluxlike["flux"] = fluxlike["optspec"] * 1
    fluxlike["uncertainty"] = fluxlike["opterr"] * 1

    return wavelike, timelike, fluxlike


def from_eureka_S3_txt(rainbow, filename, **kwargs):
    """
    Populate a Rainbow from a eureka pipeline S3 output.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filename : str
        The path to the file.
    """

    # load the Eureka event
    wavelike, timelike, fluxlike = eureadka_txt(filename)

    # populate the rainbow
    rainbow._initialize_from_dictionaries(
        wavelike=wavelike, timelike=timelike, fluxlike=fluxlike
    )
