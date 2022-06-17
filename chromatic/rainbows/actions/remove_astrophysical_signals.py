from ...imports import *


def remove_astrophysical_signals(self, method = "gradient"):
    """
    Remove astrophysical signal by creating new "removed_astrophysical_signal" fluxlike quantity using the "method" argument.

    Parameters
    ----------
    method : str
        The method used in the function.
        `gradient` will use the np.gradient to filter out long term signal (fast).

    Returns
    -------
    removed_astrophysical_signals : Rainbow
        The scatter Rainbow.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("normalize", locals())

    # TODO, think about more careful treatment of uncertainties + good/bad data
    new = self._create_copy()

    if method == "gradient":
        new.fluxlike["removed_astrophysical_signals"] = np.sqrt(2)*np.gradient(new.flux,axis=0)

    # (ignore nan warnings)
#    with warnings.catch_warnings():
#        warnings.simplefilter("ignore")
#
#        if axis.lower()[0] == "w":
#            normalization = np.nanpercentile(new.flux, percentile, axis=self.timeaxis)
#            new.fluxlike["flux"] = new.flux / normalization[:, np.newaxis]
#            try:
#                new.fluxlike["uncertainty"] = (
#                    self.uncertainty / normalization[:, np.newaxis]
#                )
#            except ValueError:
#                pass
#        elif axis.lower()[0] == "t":
#            normalization = np.nanpercentile(self.flux, percentile, axis=self.waveaxis)
#            new.fluxlike["flux"] = new.flux / normalization[np.newaxis, :]
#            try:
#                new.fluxlike["uncertainty"] = (
#                    self.uncertainty / normalization[np.newaxis, :]
#                )
#            except ValueError:
#                pass

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new
