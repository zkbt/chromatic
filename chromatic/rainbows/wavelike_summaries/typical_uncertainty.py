from ...imports import *

__all__ = ["get_typical_uncertainty"]


def get_typical_uncertainty(self, function=np.nanmedian, *args, **kwargs):
    """
    Get the typical per-wavelength uncertainty.

    Parameters
    ----------
    function : function
        What function should be used to choose the "typical"
        value for each wavelength? Good options are probably
        things like `np.nanmedian`, `np.median`, `np.nanmean`
        `np.mean`, `np.percentile`
    args : list
        Addition arguments will be passed to `function`
    kw : dict
        Additional keyword arguments will be passed to `function`

    Returns
    -------
    uncertainty_per_wavelength : np.array (wavelike)
        The uncertainty associated with each wavelength.
    """
    uncertainty_per_wavelength = function(
        self.uncertainty, *args, axis=self.timeaxis, **kwargs
    )
    return uncertainty_per_wavelength
