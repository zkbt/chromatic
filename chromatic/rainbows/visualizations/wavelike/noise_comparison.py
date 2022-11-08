from ....imports import *

__all__ = ["plot_noise_comparison"]


def plot_noise_comparison(
    self,
    ax=None,
    method="standard-deviation",
    minimum_acceptable_ok=1e-10,
    w_unit="micron",
    cmap=None,
    vmin=None,
    vmax=None,
    legend=False,
    expected_color="k",
    measured_color="auto",
    expected_label="Expected",
    measured_label="Measured",
    scatterkw={},
    plotkw={},
    filename=None,
    **kw,
):
    """
    Plot measureded per-wavelength scatter, compared to expected.

    Parameters
    ----------
    ax : Axes
        The axes into which to make this plot.
    method : string
        What method to use to obtain measured scatter. Current options are 'MAD', 'standard-deviation'.
    minimum_acceptable_ok : float
        The smallest value of `ok` that will still be included.
        (1 for perfect data, 1e-10 for everything but terrible data, 0 for all data)
    w_unit : str, Unit
        The unit for plotting wavelengths.
    cmap : str, Colormap
        The color map to use for expressing wavelength.
    vmin : Quantity
        The wavelength at the bottom of the cmap.
    vmax : Quantity
        The wavelength at the top of the cmap.
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
    legendkw : dict
        A dictionary of keywords passed to `plt.legend`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[frameon, bbox_to_anchor, ...]`
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
    **kw : dict
        Any additional keywords will be stored as `kw`.
        Nothing will happen with them.
    """

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
    plot_measured_scatter = self.get_measured_scatter(
        method=method, minimum_acceptable_ok=minimum_acceptable_ok
    )
    plot_expected_uncertainty = self.get_expected_uncertainty()
    yplots = [plot_expected_uncertainty, plot_measured_scatter]

    # define colors
    if measured_color == "auto":
        c = plot_x
        cmap = self.cmap
        norm = self.norm
        default_line_measured_color = "gray"
    else:
        c = measured_color
        default_line_measured_color = c

    plot_colors = [expected_color, default_line_measured_color]
    scatter_colors = [expected_color, c]

    labels = [expected_label, measured_label]

    for i, yplot in enumerate(yplots):

        # set default for background line plot
        this_plotkw = dict(color=plot_colors[i], label=labels[i], zorder=-1)
        this_plotkw.update(**plotkw)

        # set default for scatter plot with points
        this_scatterkw = dict(
            c=scatter_colors[i],
            cmap=cmap,
            marker="o",
            linestyle="-",
        )
        this_scatterkw.update(**scatterkw)

        # set default for error bar lines

        plt.plot(plot_x, yplot, **this_plotkw)
        plt.scatter(plot_x, yplot, **this_scatterkw)

    # add text labels to the plot
    plt.xlabel(f"Wavelength ({w_unit.to_string('latex_inline')})")
    plt.ylabel(rf"$\sigma$ ('{method}')")
    if legend:
        this_legendkw = dict(frameon=False)
        this_legendkw.update(**legendkw)
        plt.legend(**legendkw)
    plt.title(self.get("title"))

    if filename is not None:
        self.savefig(filename)
    return ax
