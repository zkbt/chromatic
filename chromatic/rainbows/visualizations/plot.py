from ...imports import *

__all__ = ["plot"]


def plot(
    self,
    ax=None,
    spacing=None,
    w_unit="micron",
    t_unit="hour",
    cmap=None,
    vmin=None,
    vmax=None,
    plotkw={},
    textkw={},
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
    plowkw : dict
        A dictionary of keywords passed to `plt.plot`
    textkw : dict
        A dictionary of keywords passed to `plt.text`
    """

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

            # grab the light curve for this particular wavelength
            lc = self.flux[i, :]

            if np.any(np.isfinite(lc)):

                # add an offset to this light curve
                plot_flux = -i * spacing + lc

                # get the color for this light curve
                color = self.get_wavelength_color(w)

                # plot the data points (with offsets)
                this_plotkw = dict(marker="o", linestyle="-", color=color)
                this_plotkw.update(**plotkw)
                plt.plot(self.time.to(t_unit), plot_flux, **this_plotkw)

                # add text labels next to each light curve
                this_textkw = dict(va="bottom", color=color)
                this_textkw.update(**textkw)
                plt.annotate(
                    f"{w.to(w_unit).value:.2f} {w_unit.to_string('latex_inline')}",
                    (min_time, 1 - (i + 0.5) * spacing),
                    **this_textkw,
                )

        # add text labels to the plot
        plt.xlabel(f"Time ({t_unit.to_string('latex_inline')})")
        plt.ylabel("Relative Flux (+ offsets)")
