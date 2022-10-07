from ....imports import *

__all__ = ["get_measured_scatter"]


def get_measured_scatter(
    self, quantity="flux", method="standard-deviation", minimum_acceptable_ok=1e-10
):
    """
    Get measured scatter for each wavelength.

    Calculate the standard deviation (or outlier-robust
    equivalent) for each wavelength, which can be compared
    to the expected per-wavelength uncertainty.

    Parameters
    ----------
    quantity : string, optional
        The `fluxlike` quantity for which we should calculate the scatter.
    method : string, optional
        What method to use to obtain measured scatter.
        Current options are 'MAD', 'standard-deviation'.
    minimum_acceptable_ok : float, optional
        The smallest value of `ok` that will still be included.
        (1 for perfect data, 1e-10 for everything but terrible data, 0 for all data)

    Returns
    -------
    scatter : array
        Wavelike array of measured scatters.
    """

    if method not in ["standard-deviation", "MAD"]:
        cheerfully_suggest(
            f"""
        '{method}' is not an available method.
        Please choose from ['MAD', 'standard-deviation'].
        """
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        scatters = np.zeros(self.nwave)
        for i in range(self.nwave):
            x, y, sigma = self.get_ok_data_for_wavelength(
                i, y=quantity, minimum_acceptable_ok=minimum_acceptable_ok
            )
            if u.Quantity(y).unit == u.Unit(""):
                y_value, y_unit = y, 1
            else:
                y_value, y_unit = y.value, y.unit
            if method == "standard-deviation":
                scatters[i] = np.nanstd(y_value)
            elif method == "MAD":
                scatters[i] = mad_std(y_value, ignore_nan=True)
        return scatters * y_unit
