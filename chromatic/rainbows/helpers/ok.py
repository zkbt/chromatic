from ...imports import *

__all__ = ["get_ok_data_for_wavelength", "get_ok_data_for_time"]


def get_ok_data_for_wavelength(
    self,
    i,
    quantity="flux",
    minimum_acceptable_ok=1,
    express_badness_with_uncertainty=False,
):
    """
    A small wrapper to get the good data from a wavelength,
    either trimming out data that are not OK or inflating the
    uncertainties to infinity.

    Parameters
    ----------
    i : int
        The wavelength index to retrieve.
    quantity : string
        Which fluxlike quantity should be retrieved? (default = 'flux')
    minimum_acceptable_ok : float
        The smallest value of `ok` that will still be included.
        (1 for perfect data, 1e-10 for everything but terrible data, 0 for all data)
    express_badness_with_uncertainty : bool
        If False, data that don't pass the `ok` cut will be removed.
        If True, data that don't pass the `ok` cut will have their
        uncertainties inflated to infinity (np.inf).

    Returns
    -------
    x : np.array
        The time.
    y : np.array
        The desired quantity (default is `flux`)
    sigma : np.array
        The uncertainty on the desired quantity
    """

    # get 1D independent variable
    x = self.time

    # get 1D array of what to keep
    ok = self.ok[i, :] >= minimum_acceptable_ok

    # get 1D array of the quantity
    y = self.get(quantity)[i, :] * 1

    # get 1D array of uncertainty
    sigma = self.fluxlike["uncertainty"][i, :] * 1

    if express_badness_with_uncertainty:
        sigma[ok == False] = np.inf
        return x, y, sigma
    else:
        return x[ok], y[ok], sigma[ok]


def get_ok_data_for_time(
    self,
    i,
    quantity="flux",
    minimum_acceptable_ok=1,
    express_badness_with_uncertainty=False,
):
    """
    A small wrapper to get the good data from a time,
    either trimming out data that are not OK or inflating the
    uncertainties to infinity.

    Parameters
    ----------
    i : int
        The time index to retrieve.
    quantity : string
        Which fluxlike quantity should be retrieved? (default = 'flux')
    minimum_acceptable_ok : float
        The smallest value of `ok` that will still be included.
        (1 for perfect data, 1e-10 for everything but terrible data, 0 for all data)
    express_badness_with_uncertainty : bool
        If False, data that don't pass the `ok` cut will be removed.
        If True, data that don't pass the `ok` cut will have their
        uncertainties inflated to infinity (np.inf).

    Returns
    -------
    x : np.array
        The time.
    y : np.array
        The desired quantity (default is `flux`)
    sigma : np.array
        The uncertainty on the desired quantity
    """

    # get 1D independent variable
    x = self.wavelength

    # get 1D array of what to keep
    ok = self.ok[:, i] >= minimum_acceptable_ok

    # get 1D array of the quantity
    y = self.get(quantity)[:, i] * 1

    # get 1D array of uncertainty
    sigma = self.fluxlike["uncertainty"][:, i] * 1

    if express_badness_with_uncertainty:
        sigma[ok == False] = np.inf
        return x, y, sigma
    else:
        return x[ok], y[ok], sigma[ok]
