from ...imports import *

__all__ = ["inject_outliers"]


def inject_outliers(self, fraction=0.01, amplitude=10):
    """
    Inject some random outliers.

    To approximately simulate cosmic rays or other
    rare weird outliers, this randomly injects
    outliers into a small fraction of pixels. For
    this simple method, outliers will have the same
    amplitude, either as a ratio above the per-data-point
    or as a fixed number (if no uncertainties exist).

    Parameters
    ----------
    fraction : float, optional
        The fraction of pixels that should get outliers.
        (default = 0.01)
    amplitude : float, optional
        If uncertainty > 0, how many sigma should outliers be?
        If uncertainty = 0, what number should be injected?
        (default = 10)

    Returns
    -------
    rainbow : Rainbow
        A new `Rainbow` object with outliers injected.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inject_outliers", locals())

    # create a copy of the existing Rainbow
    new = self._create_copy()

    # pick some random pixels to inject outliers
    outliers = np.random.uniform(0, 1, self.shape) < fraction

    # inject outliers based on uncertainty if possible
    if np.any(self.uncertainty > 0):
        new.fluxlike["injected_outliers"] = outliers * amplitude * self.uncertainty
    else:
        new.fluxlike["injected_outliers"] = outliers * amplitude

    # modify the flux
    new.fluxlike["flux"] += new.fluxlike["injected_outliers"]

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new object
    return new
