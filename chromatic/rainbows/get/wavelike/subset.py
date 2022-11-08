from ....imports import *

__all__ = [
    "get_for_wavelength",
    "get_ok_data_for_wavelength",
]


def get_for_wavelength(self, i, quantity="flux"):
    """
    Get `'quantity'` associated with wavelength `'i'`.

    Parameters
    ----------
    i : int
        The wavelength index to retrieve.
    quantity : string
        The quantity to retrieve. If it is flux-like,
        row 'i' will be returned. If it is time-like,
        the array itself will be returned.

    Returns
    -------
    quantity : array, Quantity
        The 1D array of 'quantity' corresponding to wavelength 'i'.
    """
    z = self.get(quantity)
    if np.shape(z) == self.shape:
        return z[i, :]
    elif len(z) == self.ntime:
        return z
    else:
        raise RuntimeError(
            f"""
        You tried to retrieve wavelength {i} from '{quantity}',
        but this quantity is neither flux-like nor time-like.
        It's not possible to return a time-like array. Sorry!
        """
        )


def get_ok_data_for_wavelength(
    self,
    i,
    x="time",
    y="flux",
    sigma="uncertainty",
    minimum_acceptable_ok=1,
    express_badness_with_uncertainty=False,
):
    """
    A small wrapper to get the good data from a wavelength.

    Extract a slice of data, marking data that are not `ok` either
    by trimming them out entirely or by inflating their
    uncertainties to infinity.

    Parameters
    ----------
    i : int
        The wavelength index to retrieve.
    x : string, optional
        What quantity should be retrieved as 'x'? (default = 'time')
    y : string, optional
        What quantity should be retrieved as 'y'? (default = 'flux')
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
    x : array
        The time.
    y : array
        The desired quantity (default is `flux`)
    sigma : array
        The uncertainty on the desired quantity
    """

    # get 1D independent variable
    x_values = self.get_for_wavelength(i, x) * 1

    # get 1D array of what to keep
    ok = self.ok[i, :] >= minimum_acceptable_ok

    # get 1D array of the quantity
    y_values = self.get_for_wavelength(i, y) * 1

    # get 1D array of uncertainty
    sigma_values = self.get_for_wavelength(i, sigma) * 1

    if express_badness_with_uncertainty:
        sigma_values[ok == False] = np.inf
        return x_values, y_values, sigma_values
    else:
        return x_values[ok], y_values[ok], sigma_values[ok]
