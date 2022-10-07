from ....imports import *

__all__ = ["get_median_spectrum"]


def get_median_spectrum(self):
    """
    Return a spectrum of the star, medianed over all times.

    Returns
    -------
    median_spectrum : array
        Wavelike array of fluxes.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmedian(self.flux, axis=1)
