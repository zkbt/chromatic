from ...imports import *

__all__ = [
    "setup_wavelength_colors",
    "_make_sure_cmap_is_defined",
    "get_wavelength_color",
]


def setup_wavelength_colors(self, cmap=None, vmin=None, vmax=None, log=None):
    """
    Set up a color map and normalization function for
    colors datapoints by their wavelengths.

    Parameters
    ----------
    cmap : str, Colormap
        The color map to use.
    vmin : Quantity
        The wavelength at the bottom of the cmap.
    vmax : Quantity
        The wavelength at the top of the cmap.
    log : bool
        If True, colors will scale with log(wavelength).
        If False, colors will scale with wavelength.
        If None, the scale will be guessed from the internal wscale.
    """

    # populate the cmap object
    self.cmap = plt.cm.get_cmap(cmap)

    vmin = vmin
    if vmin is None:
        vmin = np.nanmin(self.wavelength)
    vmax = vmax
    if vmax is None:
        vmax = np.nanmax(self.wavelength)

    if (self.wscale in ["log"]) or (log == True):
        self.norm = col.LogNorm(
            vmin=vmin.to("micron").value, vmax=vmax.to("micron").value
        )
    if (self.wscale in ["?", "linear"]) or (log == False):
        self.norm = col.Normalize(
            vmin=vmin.to("micron").value, vmax=vmax.to("micron").value
        )


def _make_sure_cmap_is_defined(self, cmap=None, vmin=None, vmax=None):
    """
    A helper function that can be called at the start of
    any plot that that's using wavelength-colors to make
    sure that the wavelength-based colormap has been
    defined.

    Parameters
    ----------
    cmap : str, Colormap
        The color map to use for expressing wavelength.
    vmin : Quantity
        The minimum value to use for the wavelength colormap.
    vmax : Quantity
        The maximum value to use for the wavelength colormap.
    """

    if hasattr(self, "cmap"):
        if (cmap is not None) or (vmin is not None) or (vmax is not None):
            if (
                (cmap != self.get("cmap"))
                and (vmin != self.get("vmin"))
                and (vmax != self.get("vmax"))
            ):
                cheerfully_suggest(
                    """
                It looks like you're trying to set up a new custom
                cmap and/or wavelength normalization scheme. You
                should be aware that a cmap has already been defined
                for this object; if you're visualizing the same
                rainbow in different ways, we strongly suggest
                that you not change the cmap or normalization
                between them, for visual consistency.
                """
                )
        else:
            return
    self.setup_wavelength_colors(cmap=cmap, vmin=vmin, vmax=vmax)


def get_wavelength_color(self, wavelength):
    """
    Determine the color corresponding to one or more wavelengths.

    Parameters
    ----------
    wavelength : Quantity
        The wavelength value(s), either an individual
        wavelength or an array of N wavelengths.

    Returns
    -------
    colors : array
        An array of RGBA colors [or an (N,4) array].
    """
    w_unitless = wavelength.to("micron").value
    normalized_w = self.norm(w_unitless)
    return self.cmap(normalized_w)
