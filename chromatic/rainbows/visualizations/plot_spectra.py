from ...imports import *

__all__ = ["plot_spectra"]


def plot_spectra(
    self,
    ax=None,
    spacing=0.5,
    w_unit="micron",
    t_unit="hour",
    cmap=None,
    vmin=None,
    vmax=None,
    plotkw={},
    textkw={},
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
    plotkw : dict
        A dictionary of keywords passed to `plt.plot`
    textkw : dict
        A dictionary of keywords passed to `plt.text`
    """

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
            spectrum = self.normalize().flux[:, i]

            if np.any(np.isfinite(spectrum)):

                # add an offset to this spectrum
                plot_flux = -i * spacing + spectrum

                # get the color for this spectrum
                color = "darkblue"#self.get_time_color(t)

                # plot the data points (with offsets)
                this_plotkw = dict(marker="o", linestyle="-", color= color)
                this_plotkw.update(**plotkw)
                plt.plot(self.wavelength.to(w_unit), plot_flux, **this_plotkw)

                # add text labels next to each spectrum
                this_textkw = dict(va="bottom", color='k')
                this_textkw.update(**textkw)
                plt.annotate(
                    f"{t.to(t_unit).value:.2f} {t_unit.to_string('latex_inline')}",
                    (min_wave, np.median(plot_flux) - 0.5 * spacing),
                    **this_textkw,
                )

        # add text labels to the plot
        plt.xlabel(f"Wavelength ({w_unit.to_string('latex_inline')})")
        plt.ylabel("Relative Flux (+ offsets)")
