from .plot_lightcurves import *
from .plot_spectra import *

__all__ = ["plot"]


def plot(self, xaxis="time", **kw):
    """
    Plot flux either as a sequence of offset lightcurves (default)
    or as a sequence of offset spectra.

    Parameters
    ----------
    xaxis : string
        What should be plotted on the x-axis of the plot?
        'time' will plot a different light curve for each wavelength
        'wavelength' will plot a different spectrum for each timepoint
    **kw : dict
        All other keywords will be passed along to either
        `.plot_lightcurves` or `.plot_spectra` as appropriate.
        Please see the docstrings for either of those functions
        to figure out what keyword arguments you might want to
        provide here.
    """

    if xaxis.lower()[0] == "t":
        return self.plot_lightcurves(**kw)
    elif xaxis.lower()[0] == "w":
        return self.plot_spectra(**kw)
    else:
        cheerfully_suggest("Please specify either 'time' or 'wavelength' for `.plot()`")
