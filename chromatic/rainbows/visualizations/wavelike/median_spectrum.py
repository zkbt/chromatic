from ..utilities import *
from ....imports import *

__all__ = ["plot_median_spectrum"]


def plot_median_spectrum(self, filename=None, **kw):
    """
    Plot the weighted median spectrum as a function of wavelength.

    Parameters
    ----------
    **kw : dict
        Additional keywords will be passed onward to the helper
        function `._scatter_timelike_or_wavelike`. Please see its
        docstrings for options about plot appearance and layout.
    """
    y = self.get_median_spectrum()
    self._scatter_timelike_or_wavelike(
        x=self.wavelength,
        y=y,
        ylabel=f"Flux ({_get_unit_string(y)})",
        **kw,
    )
    if filename is not None:
        self.savefig(filename)
