from ....imports import *

__all__ = ["get_expected_uncertainty"]


def get_expected_uncertainty(self, function=np.nanmedian, *args, **kw):
    """
    Get the typical per-wavelength uncertainty.

    Parameters
    ----------
    function : function, optional
        What function should be used to choose the "typical"
        value for each wavelength? Good options are probably
        things like `np.nanmedian`, `np.median`, `np.nanmean`
        `np.mean`, `np.percentile`
    *args : list, optional
        Addition arguments will be passed to `function`
    **kw : dict, optional
        Additional keyword arguments will be passed to `function`

    Returns
    -------
    uncertainty_per_wavelength : array
        The uncertainty associated with each wavelength.
    """
    uncertainty_per_wavelength = function(
        self.uncertainty, *args, axis=self.timeaxis, **kw
    )
    return uncertainty_per_wavelength
