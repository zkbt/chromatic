from ...imports import *

__all__ = ["imshow"]


def imshow(
    self,
    ax=None,
    quantity="flux",
    xaxis="time",
    w_unit="micron",
    t_unit="day",
    colorbar=True,
    aspect="auto",
    mask_ok=True,
    color_ok="tomato",
    alpha_ok=0.8,
    vmin=None,
    vmax=None,
    filename=None,
    **kw,
):
    """
    Paint a 2D image of flux as a function of time and wavelength,
    using `plt.imshow` where pixels will have constant size.

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
    aspect : str, optional
        What aspect ratio should be used for the imshow?
    mask_ok : bool, optional
        Should we mark which data are not OK?
    color_ok : str, optional
        The color to be used for masking data points that are not OK.
    alpha_ok : float, optional
        The transparency to be used for masking data points that are not OK.
    **kw : dict, optional
        All other keywords will be passed on to `plt.imshow`,
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

    # make sure some wavelength and time edges are defined
    self._make_sure_wavelength_edges_are_defined()
    self._make_sure_time_edges_are_defined()

    # set up the wavelength extent
    try:
        wmin = self.wavelength_lower[0].to_value(w_unit)
        wmax = self.wavelength_upper[-1].to_value(w_unit)
    except AttributeError:
        wmin, wmax = None, None

    if (self.wscale == "linear") and (wmin is not None) and (wmax is not None):
        wlower, wupper = wmin, wmax
        wlabel = f"{self._wave_label} ({w_unit.to_string('latex_inline')})"
    elif self.wscale == "log" and (wmin is not None) and (wmax is not None):
        wlower, wupper = np.log10(wmin), np.log10(wmax)
        wlabel = (
            r"log$_{10}$" + f"[{self._wave_label}/({w_unit.to_string('latex_inline')})]"
        )
    else:
        message = f"""
        The wavelength scale for this rainbow is '{self.wscale}',
        and there are {self.nwave} wavelength centers and
        {len(self.wavelike.get('wavelength_lower', []))} wavelength edges defined.

        It's hard to imshow something with a wavelength axis
        that isn't linearly or logarithmically uniform, or doesn't
        at least have its wavelength edges defined. We're giving up
        and just using the wavelength index as the wavelength axis.

        If you want a real wavelength axis, one solution would
        be to bin your wavelengths to a more uniform grid with
        `rainbow.bin(R=...)` (for logarithmic wavelengths) or
        `rainbow.bin(dw=...)` (for linear wavelengths)
        """
        cheerfully_suggest(message)
        wlower, wupper = -0.5, self.nwave - 0.5
        wlabel = "Wavelength Index"

    # set up the time extent
    try:
        tmin = self.time_lower[0].to_value(t_unit)
        tmax = self.time_upper[-1].to_value(t_unit)
    except AttributeError:
        tmin, tmax = None, None
    if (self.tscale == "linear") and (tmin is not None) and (tmax is not None):
        tlower, tupper = tmin, tmax
        tlabel = f"{self._time_label} ({t_unit.to_string('latex_inline')})"
    elif self.tscale == "log" and (tmin is not None) and (tmax is not None):
        tlower, tupper = np.log10(tmin), np.log10(tmax)
        tlabel = (
            r"log$_{10}$" + f"[{self._time_label}/({t_unit.to_string('latex_inline')})]"
        )
    else:
        message = f"""
        The time scale for this rainbow is '{self.tscale}',
        and there are {self.ntime} time centers and
        {len(self.timelike.get('time_lower', []))} time edges defined.

        It's hard to imshow something with a time axis
        that isn't linearly or logarithmically uniform, or doesn't
        at least have its time edges defined. We're giving up
        and just using the time index as the time axis.

        If you want a real time axis, one solution would
        be to bin your times to a more uniform grid with
        `rainbow.bin(dt=...)` (for linear times).
        """
        cheerfully_suggest(message)
        tlower, tupper = -0.5, self.ntime - 0.5
        tlabel = "Time Index"

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
        self.metadata["_imshow_extent"] = [tlower, tupper, wupper, wlower]
        xlabel, ylabel = tlabel, wlabel
        z = get_2D(quantity)
        ok = get_2D("ok")
    elif xaxis.lower()[0] == "w":
        self.metadata["_imshow_extent"] = [wlower, wupper, tupper, tlower]
        xlabel, ylabel = wlabel, tlabel
        z = get_2D(quantity).T
        ok = get_2D("ok").T
    else:
        cheerfully_suggest(
            "Please specify either `xaxis='time'` or `xaxis='wavelength'` for `.plot()`"
        )

    # figure out a good shared color limits (unless already supplied)
    vmin = vmin or np.nanpercentile(u.Quantity(z.flatten()).value * 1.0, 1)
    vmax = vmax or np.nanpercentile(u.Quantity(z.flatten()).value * 1.0, 99)

    # define some default keywords
    imshow_kw = dict(interpolation="nearest", vmin=vmin, vmax=vmax)
    imshow_kw.update(**kw)
    with quantity_support():
        plt.sca(ax)

        # create an overlaying mask of which data are OK or not
        if mask_ok:
            okimshow_kw = dict(**imshow_kw)
            okimshow_kw.update(
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
            plt.imshow(
                ok,
                extent=self.metadata["_imshow_extent"],
                aspect=aspect,
                origin="upper",
                **okimshow_kw,
            )
        plt.imshow(
            z,
            extent=self.metadata["_imshow_extent"],
            aspect=aspect,
            origin="upper",
            **imshow_kw,
        )
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
    return ax
