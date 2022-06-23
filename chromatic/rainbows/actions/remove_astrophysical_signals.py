from ...imports import *


def remove_astrophysical_signals(self, method="gradient", model=None, win_length=None, polyorder=None):
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
    
    if method == "savgol":
        if win_length is None:
            print ("You need `win_length` for this `savgol` method")
        elif polyorder is None:
            print ("You need `polyorder` for this `savgol` method")
        else:
            for i in range (0,new.nwave):
                savgolfilter = savgol_filter(new.flux[i,:],window_length=win_length,polyorder=polyorder)
                new.flux[i,:] = new.flux[i,:]/savgolfilter
    
    if method == "custom":
        if model is None:
            print ("You need fluxlike model for this `custom` method")
        elif model.shape != new.flux.shape:
            print ("Your model doesn't match flux shape")
        else:
            new.flux = new.flux/model

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new
