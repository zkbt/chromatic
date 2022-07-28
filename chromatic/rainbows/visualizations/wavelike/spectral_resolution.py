from ..utilities import *
from ....imports import *

__all__ = ["plot_spectral_resolution"]


def plot_spectral_resolution(
    self,
    pixels_per_resolution_element=1,
    filename=None,
    **kw,
):
    """
    Plot the spectral resolution as a function of wavelength.

    Parameters
    ----------
    pixels_per_resolution_element : float
        How many pixels do we consider as a resolution element?
    **kw : dict
        Additional keywords will be passed onward to the helper
        function `._scatter_timelike_or_wavelike`. Please see its
        docstrings for options about plot appearance and layout.
    """
    self._scatter_timelike_or_wavelike(
        x=self.wavelength,
        y=self.get_spectral_resolution(
            pixels_per_resolution_element=pixels_per_resolution_element
        ),
        ylabel=f"$R=\lambda/d\lambda$ ({pixels_per_resolution_element} pixel)",
        **kw,
    )
    plt.title(self.get("title"))
    if filename is not None:
        self.savefig(filename)
