from ...imports import *

__all__ = ["get_spectrum"]


def get_spectrum(self):
    """
    Return a spectrum of the star, averaged over all times.

    Parameters
    ----------


    Returns
    -------
    s : np.array (wavelike)
    """

    # FIXME, think about error weighting/masking
    return np.nanmedian(self.flux, axis=self.timeaxis)
