from ...imports import *

__all__ = ["from_text"]


def from_text(rainbow, filename, **kwargs):
    """
    Populate a Rainbow from a generic text file.
    This text file must at least contain the columns
    `wavelength`, `time`, `flux`, `uncertainty`

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated. (This is intended
        to be used as a class method of Rainbow or
        of a class derived from Rainbow, as a way of
        initializing an object from files.)
    filename : str
        The path to the file.
    """

    # load the data
    data = ascii.read(filename)

    # pull out some variables
    t = sorted(np.unique(data["time"]))
    w = sorted(np.unique(data["wavelength"]))

    fluxlike = {}
    fluxlike_keys = list(data.colnames)
    fluxlike_keys.remove("time")
    fluxlike_keys.remove("wavelength")
    for k in fluxlike_keys:
        fluxlike[k] = np.ones(shape=(len(w), len(t)))
    for i in tqdm(range(len(data)), leave=False):
        # this is slow, but general
        i_time = t == data[i]["time"]
        i_wavelength = w == data[i]["wavelength"]
        for k in fluxlike_keys:
            fluxlike[k][i_wavelength, i_time] = data[i][k]

    # do a slightly better guess of the uncertainty column
    if "uncertainty" not in fluxlike:
        for k in ["error", "flux_error", "sigma", "unc"]:
            try:
                fluxlike["uncertainty"] = fluxlike[k] * 1
                break
            except KeyError:
                pass

    timelike = {}
    try:
        t.unit
    except AttributeError:
        t *= u.day
    timelike["time"] = t

    wavelike = {}
    try:
        w.unit
    except AttributeError:
        w *= u.micron
    wavelike["wavelength"] = w

    # populate the rainbow
    rainbow._initialize_from_dictionaries(
        wavelike=wavelike, timelike=timelike, fluxlike=fluxlike
    )
