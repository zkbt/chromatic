from ...imports import *


def normalize(self, axis="wavelength", percentile=50):
    """
    Normalize by dividing through by the median spectrum and/or lightcurve.

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

    # TODO, think about more careful treatment of uncertainties + good/bad data
    new = self._create_copy()

    # (ignore nan warnings)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if axis.lower()[0] == "w":
            normalization = np.nanpercentile(new.flux, percentile, axis=self.timeaxis)
            new.fluxlike["flux"] = new.flux / normalization[:, np.newaxis]
            try:
                new.fluxlike["uncertainty"] = (
                    self.uncertainty / normalization[:, np.newaxis]
                )
            except ValueError:
                pass
        elif axis.lower()[0] == "t":
            normalization = np.nanpercentile(self.flux, percentile, axis=self.waveaxis)
            new.fluxlike["flux"] = new.flux / normalization[np.newaxis, :]
            try:
                new.fluxlike["uncertainty"] = (
                    self.uncertainty / normalization[np.newaxis, :]
                )
            except ValueError:
                pass

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new
