from ....imports import *

__all__ = ["get_average_lightcurve"]


def get_average_lightcurve(self):
    """
    Return a lightcurve of the star, averaged over all wavelengths.

    This uses `bin`, which is a horribly slow way of doing what is
    fundamentally a very simply array calculation, because we
    don't need to deal with partial pixels.

    Returns
    -------
    lightcurve : array
        Timelike array of fluxes.
    """
    return self.get_average_lightcurve_as_rainbow().flux[0, :]
