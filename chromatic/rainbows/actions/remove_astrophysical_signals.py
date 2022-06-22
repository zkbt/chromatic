from ...imports import *


def remove_astrophysical_signals(self, method = "gradient", model = None):
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
        new.flux = np.sqrt(2)*np.gradient(new.flux,axis=0)
#    if method == "highpass_filter"
#        scipy.signal.butter()
#    if method == "smooth boxcar"
#        1d or 3d
#       scipy.signal.convolve
#       lighkurve.LightCurve.flatten()
#    if method == "custom":
#        new.flux = new.flux/model

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new
