from ...imports import *

__all__ = ["get_typical_uncertainty"]


def get_typical_uncertainty(self, function=np.nanmedian):
    """
    Get the typical per-wavelength uncertainty.

    Parameters
    ----------
    function : function
        What function should be used to choose the "typical"
        value for each wavelength? Good options are probably
        things like `np.nanmedian`, `np.median`, `np.nanmean`
        `np.mean`

    Returns
    -------
    uncertainty_per_wavelength : np.array (wavelike)
        The uncertainty associated with each wavelength.
    """
    uncertainty_per_wavelength = function(self.uncertainty, axis=self.timeaxis)
    return uncertainty_per_wavelength
