from ...imports import *

__all__ = ["inflate_uncertainty"]


def inflate_uncertainty(
    self,
    method="MAD",
    remove_trends=True,
    remove_trends_method="median_filter",
    remove_trends_kw={},
    minimum_inflate_ratio=1.0,
):
    """
    Inflate uncertainties to match observed scatter.

    This is a quick and approximate tool for inflating
    the flux uncertainties in a `Rainbow` to match the
    observed scatter. With defaults, this will estimate
    the scatter using a robust median-absolute-deviation
    estimate of the standard deviation (`method='MAD'`),
    applied to time-series from which temporal trends
    have been removed (`remove_trends=True`), and inflate
    the uncertainties on a per-wavelength basis. The trend
    removal, by default by subtracting off local medians
    (`remove_trends_method='median_filter'`), will squash
    many types of both astrophysical and systematic trends,
    so this function should be used with caution in
    applicants where precise and reliable uncertainties
    are needed.

    Parameters
    ----------
    method : string
        What method to use to obtain measured scatter.
        Current options are 'MAD', 'standard-deviation'.
    remove_trends : bool
        Should we remove trends before estimating by how
        much we need to inflate the uncertainties?
    remove_trends_method : str
        What method should be used to remove trends?
        See `.remove_trends` for options.
    remove_trends_kw : dict
        What keyword arguments should be passed to `remove_trends`?
    minimum_inflate_ratio : float, optional
        the minimum inflate_ratio that can be. We don't want people
        to deflate uncertainty unless a very specific case of unstable
        pipeline output.

    Returns
    -------
    removed : Rainbow
        The Rainbow with estimated signals removed.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inflate_uncertainty", locals())

    # create a new copy
    new = self._create_copy()

    # if desired, remove trends before estimating inflation factor
    if remove_trends:
        trend_removed = new.remove_trends(**remove_trends_kw)
    else:
        trend_removed = new

    # estimate the scatter
    measured_scatter = trend_removed.get_measured_scatter(
        method=method, minimum_acceptable_ok=1e-10
    )

    # get the expected uncertainty
    expected_uncertainty = trend_removed.get_expected_uncertainty()

    # calculate the necessary inflation ratio
    inflate_ratio = measured_scatter / expected_uncertainty

    # warn if there are some inflation ratios below minimum (usually = 1)
    if np.min(inflate_ratio) < minimum_inflate_ratio:
        cheerfully_suggest(
            f"""
        {np.sum(inflate_ratio < minimum_inflate_ratio)} uncertainty inflation ratios would be below
        the `minimum_inflate_ratio` of {minimum_inflate_ratio}, so they have not been changed.
        """
        )
        inflate_ratio = np.maximum(inflate_ratio, minimum_inflate_ratio)

    # store the inflation ratio
    new.wavelike["inflate_ratio"] = inflate_ratio

    # inflate the uncertainties
    new.uncertainty = new.uncertainty * inflate_ratio[:, np.newaxis]

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new
