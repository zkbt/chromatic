from ...imports import *
from scipy.special import erfc

__all__ = ["flag_outliers"]


def flag_outliers(self, how_many_sigma=5, remove_trends=True, inflate_uncertainty=True):
    """
    Flag outliers as not `ok`.

    This examines the flux array, identifies significant outliers,
    and marks them 0 in the `ok` array. The default procedure is to use
    a median filter to remove temporal trends (`remove_trends`),
    inflate the uncertainties based on the median-absolute-deviation
    scatter (`inflate_uncertainty`), and call points outliers if they
    deviate by more than a certain number of sigma (`how_many_sigma`)
    from the median-filtered level.

    The returned `Rainbow` object should be identical to the input
    one, except for the possibility that some elements in `ok` array
    will have been marked as zero. (The filtering or inflation are
    not applied to the returned object.)

    Parameters
    ----------
    how_many_sigma : float, optional
        Standard deviations (sigmas) allowed for individual data
        points before they are flagged as outliers.
    remove_trends : bool, optional
        Should we remove trends from the flux data before
        trying to look for outliers?
    inflate_uncertainty : bool, optional
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

    # create a trend-filtered object
    if remove_trends:
        filtered = new.remove_trends(method="median_filter", size=(3, 5))
    else:
        filtered = new._create_copy()

    # update the uncertainties, if need be
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
            inflated = filtered.inflate_uncertainty(method="MAD", remove_trends=True)
    else:
        inflated = filtered

    # decide which points are outliers
    is_outlier = np.abs(inflated.flux - 1) > how_many_sigma * inflated.uncertainty

    # update the output object
    new.fluxlike["flagged_as_outlier"] = is_outlier
    new.ok = new.ok * (is_outlier == False)

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    return new
