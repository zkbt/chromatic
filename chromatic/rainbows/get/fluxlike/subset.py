from ....imports import *

__all__ = [
    "get_ok_data",
]


def get_ok_data(
    self,
    y="flux",
    minimum_acceptable_ok=1,
):
    """
    A small wrapper to get the good data as a 2D array.

    Extract fluxes as 2D array, marking data that are not `ok`
    either as nan or by inflating uncertainties to infinity.

    Parameters
    ----------
    y : array
        The desired quantity (default is `flux`)
    minimum_acceptable_ok : float, optional
        The smallest value of `ok` that will still be included.
        (1 for perfect data, 1e-10 for everything but terrible data, 0 for all data)

    Returns
    -------
    flux : array
        The fluxes, but with not-OK data replaced with nan.
    """
    # get 1D array of what to keep
    ok = self.ok >= minimum_acceptable_ok

    # create copy of flux array
    y = self.get(y) * 1

    # set bad values to nan
    y[ok == False] = np.nan

    return y
