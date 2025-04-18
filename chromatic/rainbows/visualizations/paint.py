from ...imports import *

__all__ = ["paint"]


def paint(
    self,
    ax=None,
    quantity="flux",
    xaxis="time",
    w_unit="micron",
    t_unit="day",
    colorbar=True,
    mask_ok=True,
    color_ok="tomato",
    alpha_ok=0.8,
    vmin=None,
    vmax=None,
    filename=None,
    **kw,
):
    """
    Paint a 2D image of flux as a function of time and wavelength.

    This is a wrapper that tries to choose the best between `.imshow()`
    and `.pcolormesh()` for painting a 2D map. `imshow` is faster and
    does a better job handling antialiasing for large datasets, but it
    works best for linearly or logarithmically uniform axes. `pcolormesh`
    is more flexible about non-uniform axes but is slower and doesn't
    do anything to combat antialiasing.

    Parameters
    ----------
    ax : Axes, optional
        The axes into which to make this plot.
    quantity : str, optional
        The fluxlike quantity to imshow.
        (Must be a key of `rainbow.fluxlike`).
    xaxis : str
        What to use as the horizontal axis, 'time' or 'wavelength'.
    w_unit : str, Unit, optional
        The unit for plotting wavelengths.
    t_unit : str, Unit, optional
        The unit for plotting times.
    colorbar : bool, optional
        Should we include a colorbar?
    mask_ok : bool, optional
        Should we mark which data are not OK?
    color_ok : str, optional
        The color to be used for masking data points that are not OK.
    alpha_ok : float, optional
        The transparency to be used for masking data points that are not OK.
    **kw : dict, optional
        All other keywords will be passed on to `plt.imshow` or `plt.pcolormesh`,
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[cmap, norm, interpolation, alpha, vmin, vmax]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
    """
    # record all keyword inputs in a dictionary
    inputs = locals()
    inputs.pop("self")
    kw = inputs.pop("kw")
    inputs.update(**kw)

    imshow_scales = ["linear", "log"]
    if (self.wscale in imshow_scales) and (self.tscale in imshow_scales):
        ax = self.imshow(**inputs)
    else:
        ax = self.pcolormesh(**inputs)

    return ax
