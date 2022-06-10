from ...imports import *

__all__ = ["plot_spectra"]


def plot_spectra(
    self,
    quantity="flux",
    ax=None,
    spacing=0.1,
    w_unit="micron",
    t_unit="hour",
    cmap=None,
    vmin=None,
    vmax=None,
    errorbar=False,
    scatterkw={},
    errorbarkw={},
    plotkw={},
    textkw={},
    **kw,
):
    """
    Plot flux as sequence of offset spectrum.

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
    scatterkw : dict
        A dictionary of keywords passed to `plt.scatter`
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, color, s, m, edgecolor, facecolor]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
    errorbarkw : dict
        A dictionary of keywords passed to `plt.errorbar`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[alpha, elinewidth, color, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
    plotkw : dict
        A dictionary of keywords passed to `plt.plot`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[alpha, clip_on, zorder, marker, markersize,
          linewidth, linestyle, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
    textkw : dict
        A dictionary of keywords passed to `plt.text`
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    **kw : dict
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

    min_wave = np.nanmin(self.wavelength)

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

        #  loop through times
        for i, t in enumerate(self.time):
            # grab the spectrum for this particular time
            quan = u.Quantity(self.fluxlike[quantity][:, i]).value
            fluxerr = self.fluxlike["uncertainty"][:, i]
            if np.any(np.isfinite(quan)):

                # add an offset to this spectrum
                plot_quan = -i * spacing + quan

                default_color = "black"

                # set default for background line plot
                this_plotkw = dict(color=default_color, zorder=-1)
                this_plotkw.update(**plotkw)

                # set default for scatter plot with points
                this_scatterkw = dict(
                    marker="o",
                    linestyle="-",
                    c=self.wavelength.to_value(w_unit),
                    cmap=self.cmap,
                    norm=self.norm,
                )
                this_scatterkw.update(**scatterkw)

                # set default for error bar lines
                this_errorbarkw = dict(
                    color=default_color, linewidth=0, elinewidth=1, zorder=-1
                )
                this_errorbarkw.update(**errorbarkw)

                if errorbar:
                    plt.errorbar(
                        self.wavelength.to(w_unit),
                        plot_quan,
                        yerr=fluxerr,
                        **this_errorbarkw,
                    )
                plt.plot(self.wavelength.to(w_unit), plot_quan, **this_plotkw)
                plt.scatter(self.wavelength.to(w_unit), plot_quan, **this_scatterkw)

                # add text labels next to each spectrum
                this_textkw = dict(va="bottom", color=default_color)
                this_textkw.update(**textkw)
                plt.annotate(
                    f"{t.to(t_unit).value:.2f} {t_unit.to_string('latex_inline')}",
                    (min_wave, np.median(plot_quan) - 0.5 * spacing),
                    **this_textkw,
                )

        # add text labels to the plot
        plt.xlabel(f"Wavelength ({w_unit.to_string('latex_inline')})")
        plt.ylabel("Relative Flux (+ offsets)")
        if self.get("wscale") == "log":
            plt.xscale("log")
