from ....imports import *

__all__ = ["get_average_lightcurve"]


def get_average_lightcurve(self):
    """
    Return a lightcurve of the star, averaged over all wavelengths.

    Parameters
    ----------

    Returns
    -------
    lightcurve : np.array (wavelike)
    """
    return self.get_average_lightcurve_as_rainbow().flux[0, :]
