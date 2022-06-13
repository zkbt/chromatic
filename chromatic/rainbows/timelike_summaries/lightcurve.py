from ...imports import *

__all__ = ["get_lightcurve"]


def get_lightcurve(self):
    """
    Return a lightcurve of the star, averaged over all wavelengths.

    Parameters
    ----------

    Returns
    -------
    lightcurve : np.array (wavelike)
    """
    return self.get_lightcurve_as_rainbow().flux[0, :]
