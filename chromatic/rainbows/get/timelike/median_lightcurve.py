from ....imports import *

__all__ = ["get_median_lightcurve"]


def get_median_lightcurve(self):
    """
    Return a lightcurve of the star, medianed over all wavelengths.

    Returns
    -------
    median_lightcurve : array
        Timelike array of fluxes.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmedian(self.flux, axis=0)
