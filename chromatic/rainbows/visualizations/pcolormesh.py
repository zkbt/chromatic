from ...imports import *
from ...tools.resampling import leftright_to_edges

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

    By using `.pcolormesh()`, pixels can transform based on their edges,
    so non-uniform axes are allowed. This is a tiny bit slower than
    `.imshow()`, but otherwise very similar. `.paint()` will try to
    choose the best between `.imshow()` and `.pcolormesh()`.

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
    vmin = vmin or np.nanpercentile(remove_unit(z).flatten(), 1)
    vmax = vmax or np.nanpercentile(remove_unit(z).flatten(), 99)

    # define some default keywords
    pcolormesh_kw = dict(shading="flat", vmin=vmin, vmax=vmax)
    # I want to include `antialiased=True` to render nicely whether the plotting
    # pixels super- or sub-sample the data pixels, but it appears not to have
    # an effect in `pcolormesh`, so we just warn the user instead.

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
                remove_unit(x),
                remove_unit(y),
                remove_unit(ok),
                **okpcolormesh_kw,
            )
        pcolormeshed = plt.pcolormesh(
            remove_unit(x),
            remove_unit(y),
            remove_unit(z),
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

        # offer warning if aliasing might be a problem
        plt.draw()
        bbox = pcolormeshed.get_window_extent()
        render_pixels = np.abs([bbox.width, bbox.height])
        data_pixels = np.array(z.shape[::-1])
        ratios = render_pixels / data_pixels
        aliasing_warning_threshold = 2
        if np.min(ratios) < aliasing_warning_threshold:
            cheerfully_suggest(
                f"""
            In using`.pcolormesh`, [display pixels] / [data pixels] =
            {np.round(render_pixels)} / {data_pixels} = {ratios} {(xlabel.split()[0], ylabel.split()[0])}
            Is less than the suggested threshold of {aliasing_warning_threshold}.

            This suggests that aliasing/moiré might be a problem, where too many
            data pixels are trying to be displayed with too few pixels, and the
            choices `matplotlib` makes for how to do that might not be intuitive.

            Here are possible solutions:
                - Use `.bin()` to decrease the number of data pixels in wavelength
                and/or time, effectively averaging before displaying, rather than
                asking `matplotlib` to decide how to visually average adjacent data.
                - Increase the `dpi` of the figure in which this appears, so there are
                enough display pixels to represent all the data pixels being shown.
                - Use `.imshow()` instead of `.pcolormesh`, which can do better
                built-in handling of antialiasing for large data arrays. Since `.imshow()`
                can only show uniform wavelength and time grids, non-uniform grids will
                be labeled via index instead of actual value.
            """
            )
    if filename is not None:
        self.savefig(filename)
    return ax
