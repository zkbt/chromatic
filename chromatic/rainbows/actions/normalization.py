from ...imports import *

__all__ = ["normalize", "_is_probably_normalized"]


def normalize(self, axis="wavelength", percentile=50):
    """
    Normalize by dividing through by the median spectrum and/or lightcurve.

    This normalizes a `Rainbow` by estimating dividing
    through by a wavelength-dependent normalization. With
    default inputs, this would normalize each wavelength
    to have flux values near 1, to make it easier to see
    differences across time (such as a transit or eclipse).
    This function could also be used to divide through by
    a median light curve, to make it easier to see variations
    across wavelength.

    Parameters
    ----------
    axis : str
        The axis that should be normalized out.
        `w` or `wave` or `wavelength` will divide out the typical spectrum.
        `t` or `time` will divide out the typical light curve

    percentile : float
        A number between 0 and 100, specifying the percentile
        of the data along an axis to use as the reference.
        The default of `percentile=50` corresponds to the median.
        If you want to normalize to out-of-transit, maybe you
        want a higher percentile. If you want to normalize to
        the baseline below a flare, maybe you want a lower
        percentage.

    Returns
    -------
    normalized : Rainbow
        The normalized Rainbow.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("normalize", locals())

    # create an empty copy
    new = self._create_copy()

    # shortcut for the first letter of the axis
    a = axis.lower()[0]

    # (ignore nan warnings)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # get fluxes, with not-OK replaced with nans
        flux_for_normalizing = new.get_ok_data()
        negative_normalization_message = ""
        if a == "w":
            normalization = np.nanpercentile(
                flux_for_normalizing, percentile, axis=self.timeaxis
            )
            for k in self._keys_that_respond_to_math:
                new.fluxlike[k] = new.get(k) / normalization[:, np.newaxis]
            try:
                new.fluxlike["uncertainty"] = (
                    self.uncertainty / normalization[:, np.newaxis]
                )
            except ValueError:
                pass

        elif a == "t":
            normalization = np.nanpercentile(
                flux_for_normalizing, percentile, axis=self.waveaxis
            )
            for k in self._keys_that_respond_to_math:
                new.fluxlike[k] = new.get(k) / normalization[np.newaxis, :]
            try:
                new.fluxlike["uncertainty"] = (
                    self.uncertainty / normalization[np.newaxis, :]
                )
            except ValueError:
                pass

    if a in "wt":
        thing = {"w": "wavelengths", "t": "times"}[a]
        fix = {
            "w": """
                ok = rainbow.get_median_spectrum() > 0
                rainbow[ok, :].normalize()
        """,
            "t": """
                ok = rainbow.get_median_lightcurve() > 0
                rainbow[:, ok].normalize()
        """,
        }[a]
        if np.any(normalization < 0):
            cheerfully_suggest(
                f"""
            There are {np.sum(normalization < 0)} negative {thing} that
            are going into the normalization of this Rainbow. If you're
            not expecting negative fluxes, it may be useful to trim them
            away with something like:

            {fix}

            Otherwise, watch out that your fluxes and uncertainties may
            potentially have flipped sign!
            """
            )

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new


def _is_probably_normalized(
    self,
):
    """
    A helper to guess whether this `Rainbow` has been normalized or not.
    """

    # was there a normalization step?
    is_normalized = "normalize" in self.history()

    # are values generally close to 1?
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spectrum = self.get_median_spectrum()

        sigma = np.maximum(
            u.Quantity(self.get_expected_uncertainty()).value,
            u.Quantity(self.get_measured_scatter(method="MAD")).value,
        )
    try:
        sigma_value = u.Quantity(sigma).value
        spectrum_value = u.Quantity(spectrum).value
        assert np.any(sigma_value > 0)
        is_close = np.nanpercentile(np.abs(spectrum_value - 1) / sigma, 95) < 5
    except AssertionError:
        is_close = np.nanpercentile(np.abs(spectrum_value - 1), 95) < 0.1
    return is_normalized or is_close
