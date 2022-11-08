from ....imports import *

__all__ = ["plot_noise_comparison_in_bins"]


def plot_noise_comparison_in_bins(
    self,
    ax=None,
    cmap=None,
    vmin=None,
    vmax=None,
    expected=True,
    measured_errorbarkw={},
    expected_plotkw={},
    filename=None,
    ylim=[1e-5, 1e-1],
    **kw,
):
    """
    ax : Axes
        The axes into which the plot should be drawn.
    cmap : str, Colormap
        The color map to use for expressing wavelength.
    vmin : Quantity
        The wavelength at the bottom of the cmap.
    vmax : Quantity
        The wavelength at the top of the cmap.
    errorbar_kw : dict

    **kw : dictionary
        Additional keywords are passed to `.get_measured_scatter_in_bins`.
        Please see the docstring for that method for options.

    """

    # figure out where the plot should be drawn
    if ax is None:
        ax = plt.subplot()
    plt.sca(ax)

    # make sure the color map is set up
    self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

    # calculate the binned down scatters
    x = self.get_measured_scatter_in_bins(**kw)

    # loop through wavelengths
    for i, w in enumerate(self.wavelength):
        color = self.get_wavelength_color(w)
        ekw = dict(marker="o", color=color)
        ekw.update(**measured_errorbarkw)
        plt.errorbar(x["N"], x["scatters"][i], x["uncertainty"][i], **ekw)
        if expected:
            pkw = dict(color=color, alpha=0.2)
            pkw.update(**expected_plotkw)
            plt.plot(x["N"], x["expectation"][i], **pkw)

    # set good plot appearance default
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim(*ylim)
    plt.xlim(1, np.nanmax(x["N"]))
    plt.xlabel("# of Times in Bin")

    t_unit = u.minute
    twin_ax = ax.twiny()

    plt.sca(twin_ax)
    dt = x["dt"].to(t_unit).value
    plt.plot(dt, x["expectation"][i], alpha=0)
    plt.xlim(np.nanmin(dt), np.nanmax(dt))
    plt.xscale("log")
    plt.xlabel(f"Duration of Bin ({t_unit.to_string()})")

    plt.sca(ax)
    # plt.title(self.get("title"))

    if filename is not None:
        self.savefig(filename)
    return ax
