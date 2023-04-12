from ....imports import *

__all__ = [
    "get_ok_data",
]


def get_ok_data(
    self,
    xy="flux",
    sigma="uncertainty",
    minimum_acceptable_ok=1,
    express_badness_with_uncertainty=False,
):
    """
    A small wrapper to get the good data from flux.

    Extract a subset of data, marking data that are not `ok` either
    by trimming them out entirely or by inflating their
    uncertainties to infinity.

    Parameters
    ----------
    xy : string, optional
        What 2D quantity should be retrieved as 'xy'? (default = 'flux')
    sigma : string, optional
        What quantity should be retrieved as 'sigma'? (default = 'uncertainty')
    minimum_acceptable_ok : float, optional
        The smallest value of `ok` that will still be included.
        (1 for perfect data, 1e-10 for everything but terrible data, 0 for all data)
    express_badness_with_uncertainty : bool, optional
        If False, data that don't pass the `ok` cut will be removed.
        If True, data that don't pass the `ok` cut will have their
        uncertainties inflated to infinity (np.inf).

    Returns
    -------
    xy : array
        The desired quantity (default is `flux`)
    sigma : array
        The uncertainty on the desired quantity
    """

    # get 2D variable
    # x_values = self.get_for_wavelength(i, x) * 1
    xy_values = self.get(xy)

    # get 2D array of what to keep
    ok = self.ok >= minimum_acceptable_ok

    # get 2D array of uncertainty
    sigma_values = self.get(sigma)

    if express_badness_with_uncertainty:
        sigma_values[ok == False] = np.inf
        return xy_values, sigma_values
    else:
        warnings.warn("WARNING: You are getting a 1D version of a 2D array!")
        return xy_values[ok], sigma_values[ok]
