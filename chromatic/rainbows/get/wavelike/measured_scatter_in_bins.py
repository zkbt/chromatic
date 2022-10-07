from ....imports import *

__all__ = ["get_measured_scatter_in_bins"]


def get_measured_scatter_in_bins(
    self, ntimes=2, nbins=4, method="standard-deviation", minimum_acceptable_ok=1e-10
):
    """
    Get measured scatter in time bins of increasing sizes.

    For uncorrelated Gaussian noise, the scatter should
    decrease as 1/sqrt(N), where N is the number points
    in a bin. This function calculates the scatter for
    a range of N, thus providing a quick test for
    correlated noise.

    Parameters
    ----------
    ntimes : int
        How many times should be binned together? Binning will
        continue recursively until fewer that nbins would be left.
    nbins : int
        What's the smallest number of bins that should be used to
        calculate a scatter? The absolute minimum is 2.
    method : string
        What method to use to obtain measured scatter. Current options are 'MAD', 'standard-deviation'.
    minimum_acceptable_ok : float
        The smallest value of `ok` that will still be included.
        (1 for perfect data, 1e-10 for everything but terrible data, 0 for all data)

    Returns
    -------
    scatter_dictionary : dict
        Dictionary with lots of information about scatter in bins per wavelength.
    """

    from ...rainbow import Rainbow

    if "remove_trends" in self.history():
        cheerfully_suggest(
            f"""
        The `remove_trends` function was applied to this `Rainbow`,
        making it very plausible that some long-timescale signals
        and/or noise have been suppressed. Be suspicious of binned
        scatters on long timescales.
        """
        )

    # create a simplified rainbow so we don't waste time binning
    simple = Rainbow(
        time=self.time,
        wavelength=self.wavelength,
        flux=self.flux,
        uncertainty=self.uncertainty,
        ok=self.ok,
    )

    # loop through binning until done
    binnings = [simple]
    N = [1]
    while binnings[-1].ntime > ntimes * nbins:
        binnings.append(
            binnings[-1].bin(ntimes=ntimes, minimum_acceptable_ok=minimum_acceptable_ok)
        )
        N.append(N[-1] * ntimes)

    scatters = [b.get_measured_scatter(method=method) for b in binnings]
    expectation = [b.get_expected_uncertainty() for b in binnings]
    uncertainty_on_scatters = (
        scatters
        / np.sqrt(2 * (np.array([b.ntime for b in binnings]) - 1))[:, np.newaxis]
    )
    dt = [np.median(np.diff(b.time)) for b in binnings]

    return dict(
        N=np.array(N),
        dt=u.Quantity(dt),
        scatters=np.transpose(scatters),
        expectation=np.transpose(expectation),
        uncertainty=np.transpose(uncertainty_on_scatters),
    )
    # (see equation 3.48 of Sivia and Skilling for the uncertainty on sigma)
