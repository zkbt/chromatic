from ..utilities import *
from ....imports import *

__all__ = ["plot_average_spectrum"]


def plot_average_spectrum(self, **kw):
    """
    Plot the weighted average average_average_spectrum as a function of wavelength.

    Parameters
    ----------
    **kw : dict
        Additional keywords will be passed onward to the helper
        function `._scatter_timelike_or_wavelike`. Please see its
        docstrings for options about plot appearance and layout.
    """
    y = self.get_average_spectrum()
    y_unit = u.Quantity(y).unit
    if y_unit == u.Unit(""):
        unit_string = "unitless"
    else:
        unit_string = y_unit.to_string("latex_inline")
    self._scatter_timelike_or_wavelike(
        x=self.wavelength,
        y=y,
        ylabel=f"Flux ({unit_string})",
        **kw,
    )
