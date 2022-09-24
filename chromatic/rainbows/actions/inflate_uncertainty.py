from ...imports import *

__all__ = ["inflate_uncertainty"]


def inflate_uncertainty(
    self,
    minimum_inflate_ratio = 1.0,
    **kw,
):  
    """
    A quick tool to inflate uncertainties to photon noise.

    Parameters
    ----------
    method : str
        What method should be used to make an approximate model
        for smooth trends that will then be subtracted off?

        `differences` will do an extremely rough filtering
        of replacing the fluxes with their first differences.
        Trends that are smooth relative to the noise will
        be removed this way, but sharp features will remain.
        Required keywords:
            None.

        `median_filter` is a wrapper for scipy.signal.median_filter.
        It smoothes each data point to the median of its surrounding
        points in time and/or wavelength. Required keywords:
            `size` = centered on each point, what shape rectangle
            should be used to select surrounding points for median?
            The dimensions are (nwavelengths, ntimes), so `size=(3,7)`
            means we'll take the median across three wavelengths and
            seven times. Default is `(1,5)`.

        `savgol_filter` is a wrapper for scipy.signal.savgol_filter.
        It applies a Savitzky-Golay filter for polynomial smoothing.
        Required keywords:
            `window_length` = the length of the filter window,
            which must be a positive odd integer. Default is `5`.
            `polyorder` = the order of the polynomial to use.
            Default is `2`.

        `polyfit` is a wrapper for numpy.polyfit to use a weighted
        linear least squares polynomial fit to remove smooth trends
        in time. Required keywods:
            `deg` = the polynomial degree, which must be a positive
            integer. Default is `1`, meaning a line.

        `custom` allow users to pass any fluxlike array of model
        values for an astrophysical signal to remove it. Required
        keywords:
            `model` = the (nwavelengths, ntimes) model array

    minimum_inflate_ratio : float
	the minimum inflate_ratio that can be. We don't want people 
	to deflate uncertainty unless a very specific case of unstable 
	pipeline output.

    kw : dict
        Any additional keywords will be passed to the function
        that does the filtering. See `method` keyword for options.

    Returns
    -------
    removed : Rainbow
        The Rainbow with estimated signals removed.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inflate_uncertainty", locals())

    # TODO, think about more careful treatment of uncertainties + good/bad data
    new = self._create_copy()
	    
    trend_removed = new.remove_trends(**kw)
    measured_scatter = trend_removed.get_measured_scatter(method='MAD',minimum_acceptable_ok=1e-10)

    expected_uncertainty = trend_removed.get_expected_uncertainty()

    inflate_ratio = measured_scatter/expected_uncertainty
    new.wavelike['inflate_ratio'] = inflate_ratio
    if np.min(inflate_ratio) < minimum_inflate_ratio:
    	raise Exception(f"One or more inflate ratios are below the minimum_inflate_ratio = {minimum_inflate_ratio}, if you really want to proceed please consider lower minimum_inflate_ratio")

    new.uncertainty = new.uncertainty*inflate_ratio[:,np.newaxis]  
  
    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new
