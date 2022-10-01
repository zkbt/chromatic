from ...imports import *
from scipy.special import erfc

__all__ = ["flag_outliers"]


def flag_outliers(self, how_many_sigma=5, inflate_uncertainty=True):
    """
    Flag outliers on a chromatic rainbow.

    Parameters
    ----------
    how_many_sigma : float
        Standard deviations (sigmas) allowed for individual data
        points before they are flagged as outliers.
    inflate_uncertainty : bool
        Should uncertainties per wavelength be inflated to
        match the (MAD-based) standard deviation of the data?

    Returns
    -------
    rainbow : Rainbow
        A new Rainbow object with the outliers flagged as 0 in `.ok`
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("flag_outliers", locals())

    # create a copy of the existing rainbow
    new = self._create_copy()

    # how many outliers are expected from noise alone
    outliers_expected_from_normal_distribution = erfc(how_many_sigma) * self.nflux * 2
    if outliers_expected_from_normal_distribution >= 1:
        cheerfully_suggest(
            f"""
        When drawing from a normal distribution, an expected {outliers_expected_from_normal_distribution:.1f} out of
        the total {self.nflux} datapoints in {self} would be marked
        as a >{how_many_sigma} sigma outlier.

        If you don't want to accidentally clip legitimate data points that
        might have arisen merely by chance, please consider setting the
        outlier flagging threshold (`sigma=`) to a larger value.
        """
        )

    # create a trend-filtered object and use it to find outliers
    filtered = new.remove_trends(method="median_filter", size=(3, 5))
    if np.all(filtered.uncertainty == 0):
        filtered.uncertainty = (
            np.ones(filtered.shape)
            * filtered.get_measured_scatter(method="MAD")[:, np.newaxis]
        )
        inflate_uncertainty = False

    # inflate the per-wavelength uncertainties, as needed
    if inflate_uncertainty:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            inflated = filtered.inflate_uncertainty(method="MAD")
    else:
        inflated = filtered
    is_outlier = np.abs(inflated.flux - 1) > how_many_sigma * inflated.uncertainty

    # update the output object
    new.fluxlike["flagged_as_outlier"] = is_outlier
    new.ok = new.ok * (is_outlier == False)

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    return new
