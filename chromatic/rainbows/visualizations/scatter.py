from ...imports import *

__all__ = ["scatter"]


def scatter(
    self,
    ax=None,
    quantity="flux",
    xaxis="time",
    w_unit="micron",
    t_unit="day",
    s=9,
    colorbar=True,
    mask_ok=True,
    vmin=None,
    vmax=None,
    filename=None,
    **kw,
):
    """
    Paint a sparse 2D image of flux as a function of time and wavelength,
    using `plt.scatter` to place points colored by their flux.

    Parameters
    ----------
    ax : Axes, optional
        The axes into which to make this plot.
    quantity : str, optional
        The fluxlike quantity to scatter.
        (Must be a key of `rainbow.fluxlike`).
    w_unit : str, Unit, optional
        The unit for plotting wavelengths.
    t_unit : str, Unit, optional
        The unit for plotting times.
    colorbar : bool, optional
        Should we include a colorbar?
    mask_ok : bool, optional
        Should we hide data that are not OK?
    color_ok : str, optional
        The color to be used for masking data points that are not OK.
    alpha_ok : float, optional
        The transparency to be used for masking data points that are not OK.
    **kw : dict, optional
        All other keywords will be passed on to `plt.scatter`,
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[cmap, norm, interpolation, alpha, vmin, vmax]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
    """

    # self.speak(f'imshowing')
    if ax is None:
        ax = plt.subplot()

    # get units
    w_unit, t_unit = u.Unit(w_unit), u.Unit(t_unit)

    # set up the x and y axes
    w, t = np.meshgrid(self.wavelength, self.time, indexing="ij")
    wlabel = f"{self._wave_label} ({w_unit.to_string('latex_inline')})"
    tlabel = f"{self._time_label} ({t_unit.to_string('latex_inline')})"
    if (self.tscale == "linear" or self.tscale == "log") and (
        self.wscale == "linear" or self.wscale == "log"
    ):
        cheerfully_suggest(
            f"""
        Times are {self.tscale}ly spaced, and wavelengths are {self.wscale}ly spaced.
        It's likely that `.imshow` or `.pcolormesh` are better ways to display this Rainbow.
        """
        )

    def get_2D(k):
        """
        A small helper to get a 2D quantity. This is a bit of
        a kludge to help with weird cases of duplicate keys
        (for example where 'wavelength' might appear in both
        `wavelike` and `fluxlike`).
        """
        z = self.get(k)
        if np.shape(z) == self.shape:
            return z
        else:
            return self.fluxlike.get(k, None)

    z = get_2D(quantity)
    ok = get_2D("ok")

    # hide bad data if requested
    if mask_ok:
        w = w[ok]
        t = t[ok]
        z = z[ok]

    # choose between time and wavelength on the x-axis
    if xaxis.lower()[0] == "t":
        xlabel, ylabel = tlabel, wlabel
        x, y = t, w
    elif xaxis.lower()[0] == "w":
        xlabel, ylabel = wlabel, tlabel
        x, y = w, t
    else:
        cheerfully_suggest(
            "Please specify either `xaxis='time'` or `xaxis='wavelength'` for `.scatter()`"
        )

    # figure out a good shared color limits (unless already supplied)
    vmin = vmin or np.nanpercentile(u.Quantity(z.flatten()).value * 1.0, 1)
    vmax = vmax or np.nanpercentile(u.Quantity(z.flatten()).value * 1.0, 99)

    # define some default keywords
    scatter_kw = dict(s=s, vmin=vmin, vmax=vmax)
    scatter_kw.update(**kw)
    with quantity_support():
        plt.sca(ax)
        plt.scatter(x, y, c=z, **scatter_kw)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        if colorbar:
            plt.colorbar(
                ax=ax,
                label=u.Quantity(z).unit.to_string("latex_inline"),
            )
        plt.title(self.get("title"))

    if filename is not None:
        self.savefig(filename)
