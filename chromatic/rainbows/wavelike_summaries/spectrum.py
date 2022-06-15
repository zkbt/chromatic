from ...imports import *

__all__ = ["get_spectrum"]


def get_spectrum(self):
    """
    Return a spectrum of the star, averaged over all times.

    Parameters
    ----------

    Returns
    -------
    spectrum : np.array (wavelike)
    """
    return self.get_spectrum_as_rainbow().flux[:, 0]
