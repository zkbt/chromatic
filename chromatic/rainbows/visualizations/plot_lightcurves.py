from ...imports import *

__all__ = ["plot_lightcurves"]


def plot_lightcurves(
    self,
    quantity="flux",
    ax=None,
    spacing=None,
    w_unit="micron",
    t_unit="hour",
    cmap=None,
    vmin=None,
    vmax=None,
    errorbar=False,
    plotkw={},
    errorbarkw={},
    textkw={},
    **kw,
):
    """
    Plot flux as sequence of offset light curves.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    spacing : None, float
        The spacing between light curves.
        (Might still change how this works.)
        None uses half the standard dev of entire flux data.
    w_unit : str, astropy.unit.Unit
        The unit for plotting wavelengths.
    t_unit : str, astropy.unit.Unit
        The unit for plotting times.
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength.
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap.
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap.
    errorbar : boolean
        Should we plot errorbars?
    plotkw : dict
        A dictionary of keywords passed to `plt.plot`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[alpha, clip_on, zorder, marker, markersize,
          linewidth, linestyle, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
    errorbarkw : dict
        A dictionary of keywords passed to `plt.errorbar`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[alpha, elinewidth, color, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
    textkw : dict
        A dictionary of keywords passed to `plt.text`
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    kw : dict
        Any additional keywords will be stored as `kw`.
        Nothing will happen with them.
    """
    if len(kw) > 0:
        message = f"""
        You provided the keyword argument(s)
        {kw}
        but this function doesn't know how to
        use them. Sorry!
        """
        warnings.warn(message)

    # make sure that the wavelength-based colormap is defined
    self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

    w_unit, t_unit = u.Unit(w_unit), u.Unit(t_unit)

    min_time = np.nanmin(self.time)

    # make sure ax is set up
    if ax is None:
        ax = plt.subplot()
    plt.sca(ax)

    # figure out the spacing to use
    if spacing is None:
        try:
            spacing = ax._most_recent_chromatic_plot_spacing
        except AttributeError:
            spacing = 3 * np.nanstd(self.flux)
    ax._most_recent_chromatic_plot_spacing = spacing

    # TO-DO: check if this Rainbow has been normalized
    '''warnings.warn(
        """
    It's not clear if/how this object has been normalized.
    Be aware that the baseline flux levels may therefore
    be a little bit funny in .plot()."""
    )'''
    with quantity_support():

        #  loop through wavelengths
        for i, w in enumerate(self.wavelength):

            # grab the quantity and yerr for this particular wavelength
            quan = self.fluxlike[quantity][i, :]
            uncertainty = self.fluxlike["uncertainty"][i, :]

            if np.any(np.isfinite(quan)):

                # add an offset to this quantity
                plot_quan = -i * spacing + quan

                # get the color for this quantity
                color = self.get_wavelength_color(w)

                # plot the data points (with offsets)
                this_plotkw = dict(marker="o", linestyle="-", color=color)
                this_plotkw.update(**plotkw)

                # set default for error bar lines
                this_errorbarkw = dict(
                    color=color, linewidth=0, elinewidth=1, zorder=-1
                )
                this_errorbarkw.update(**errorbarkw)

                if errorbar:
                    plt.errorbar(
                        self.time.to(t_unit),
                        plot_quan,
                        yerr=uncertainty,
                        **this_errorbarkw,
                    )
                plt.plot(self.time.to(t_unit), plot_quan, **this_plotkw)

                # add text labels next to each quantity plot
                this_textkw = dict(va="bottom", color=color)
                this_textkw.update(**textkw)
                plt.annotate(
                    f"{w.to(w_unit).value:.2f} {w_unit.to_string('latex_inline')}",
                    (min_time, np.median(plot_quan) - 0.5 * spacing),
                    **this_textkw,
                )

        # add text labels to the plot
        plt.xlabel(f"{self._time_label} ({t_unit.to_string('latex_inline')})")
        plt.ylabel("Relative Flux (+ offsets)")
