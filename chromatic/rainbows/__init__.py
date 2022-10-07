from .withmodel import *
from .rainbow import *
from .simulated import *
from .multi import *
from .writers import *


def read_rainbow(filepath, **kw):
    """
    A friendly wrapper to load time-series spectra and/or
    multiwavelength light curves into a `chromatic` Rainbow
    object. It will try its best to pick the best reader
    and return the most useful kind of object.
    🦋🌅2️⃣🪜🎬👀🇮🇹📕🧑‍🏫🌈

    Parameters
    ----------
    filepath : str, list
        The file or files to open.
    **kw : dict
        All other keyword arguments will be passed to
        the `Rainbow` initialization.

    Returns
    -------
    rainbow : Rainbow, RainbowWithModel
        The loaded data!
    """
    r = Rainbow(filepath, **kw)
    if "model" in r.fluxlike:
        return RainbowWithModel(**r._get_core_dictionaries())
    else:
        return r
