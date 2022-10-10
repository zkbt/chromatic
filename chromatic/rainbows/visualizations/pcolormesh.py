from ...imports import *
from ...resampling import leftright_to_edges

__all__ = ["pcolormesh"]


def pcolormesh(
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

    By using `.pcolormesh`, pixels can transform based on their edges,
    so non-uniform axes are allowed. This is a tiny bit slower than
    `.imshow`, but otherwise very similar.

    Parameters
    ----------
    ax : Axes, optional
        The axes into which to make this plot.
    quantity : str, optional
        The fluxlike quantity to imshow.
        (Must be a key of `rainbow.fluxlike`).
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
        All other keywords will be passed on to `plt.pcolormesh`,
        so you can have more detailed control over the plot
        appearance. Common keyword argumentsvli might include:
        `[cmap, norm, alpha, vmin, vmax]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pcolormesh.html
    """

    # self.speak(f'imshowing')
    if ax is None:
        ax = plt.subplot()

    # get units
    w_unit, t_unit = u.Unit(w_unit), u.Unit(t_unit)

    # make sure some wavelength and time edges are defined
    self._make_sure_wavelength_edges_are_defined()
    self._make_sure_time_edges_are_defined()

    # set up the wavelength and time edges
    w_edges = leftright_to_edges(
        self.wavelength_lower.to_value(w_unit), self.wavelength_upper.to_value(w_unit)
    )
    t_edges = leftright_to_edges(
        self.time_lower.to_value(t_unit), self.time_upper.to_value(t_unit)
    )

    wlabel = f"{self._wave_label} ({w_unit.to_string('latex_inline')})"
    tlabel = f"{self._time_label} ({t_unit.to_string('latex_inline')})"

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

    if xaxis.lower()[0] == "t":
        x, y = t_edges, w_edges
        xlabel, ylabel = tlabel, wlabel
        z = get_2D(quantity)
        ok = get_2D("ok")
    elif xaxis.lower()[0] == "w":
        x, y = w_edges, t_edges
        xlabel, ylabel = wlabel, tlabel
        z = get_2D(quantity).T
        ok = get_2D("ok").T
    else:
        cheerfully_suggest(
            "Please specify either `xaxis='time'` or `xaxis='wavelength'` for `.plot()`"
        )

    # figure out a good shared color limits (unless already supplied)
    vmin = vmin or np.nanpercentile(u.Quantity(z.flatten()).value, 1)
    vmax = vmax or np.nanpercentile(u.Quantity(z.flatten()).value, 99)

    # define some default keywords
    pcolormesh_kw = dict(shading="flat", vmin=vmin, vmax=vmax)
    pcolormesh_kw.update(**kw)
    with quantity_support():
        plt.sca(ax)
        if mask_ok:
            okpcolormesh_kw = dict(**pcolormesh_kw)
            okpcolormesh_kw.update(
                cmap=one2another(
                    bottom=color_ok,
                    top=color_ok,
                    alpha_bottom=alpha_ok,
                    alpha_top=0,
                ),
                zorder=10,
                vmin=0,
                vmax=1,
            )
            plt.pcolormesh(
                x,
                y,
                ok,
                **okpcolormesh_kw,
            )
        plt.pcolormesh(
            x,
            y,
            z,
            **pcolormesh_kw,
        )
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        if colorbar:
            plt.colorbar(
                ax=ax,
                label=u.Quantity(z).unit.to_string("latex_inline"),
            )
        # emulate origin = upper for imshow (y starts at top)
        plt.ylim(y[-1], y[0])
        plt.title(self.get("title"))

    if filename is not None:
        self.savefig(filename)
    return ax
