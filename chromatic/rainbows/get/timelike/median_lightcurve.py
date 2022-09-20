from ....imports import *

__all__ = ["get_median_lightcurve"]


def get_median_lightcurve(self):
    """
    Return a lightcurve of the star, medianed over all wavelengths.

    Parameters
    ----------

    Returns
    -------
    median_lightcurve : np.array (wavelike)
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmedian(self.flux, axis=0)
