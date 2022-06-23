from ...imports import *

__all__ = ["plot_noise_comparison"]


def plot_noise_comparison(
    self,
    ax=None,
    w_unit="micron",
    cmap=None,
    norm=None,
    vmin=None,
    vmax=None,
    typical_uncertainty_color='k',
    measured_scatter_color='auto',
    typical_uncertainty_label = "Typical Uncertainty",
    measured_scatter_label = "SD of the Flux",
    scatterkw={},
    plotkw={},
    **kw,
):
    """
    Plot flux as sequence of offset spectrum.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    w_unit : str, astropy.unit.Unit
        The unit for plotting wavelengths.
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength.
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

    w_unit = u.Unit(w_unit)

    # make sure ax is set up
    if ax is None:
        ax = plt.subplot()
    plt.sca(ax)
    
    # set wavelength as x-axis for the plot
    w = self.wavelength
    plot_x = w.to_value(w_unit)

    # specify y_plots
    plot_measured_scatter = self.get_measured_scatter()
    plot_typical_uncertainty = self.get_typical_uncertainty()
    yplots = [self.get_measured_scatter(),self.get_typical_uncertainty()]

    # define colors
    if measured_scatter_color == "auto":
        c = plot_x
        cmap = self.cmap
        norm = self.norm
    else:
        c = measured_scatter_color

    default_line_measured_scatter_color = "red"

    plot_colors = [typical_uncertainty_color, default_line_measured_scatter_color]
    scatter_colors = [typical_uncertainty_color, c]

    labels = [typical_uncertainty_label,measured_scatter_label]

    for i,yplot in enumerate(yplots):

        # set default for background line plot
        this_plotkw = dict(color=plot_colors[i], label=labels[i], zorder=-1)
        this_plotkw.update(**plotkw)

        # set default for scatter plot with points
        this_scatterkw = dict(
            c = scatter_colors[i],
            cmap = cmap,
            norm = norm,
            marker="o",
            linestyle="-",
        )
        this_scatterkw.update(**scatterkw)

        # set default for error bar lines

        plt.plot(plot_x, yplot, **this_plotkw)
        plt.scatter(plot_x, yplot,**this_scatterkw)


    # add text labels to the plot
    plt.xlabel(f"Wavelength ({w_unit.to_string('latex_inline')})")
    plt.ylabel("Relative Flux")
    plt.legend()
