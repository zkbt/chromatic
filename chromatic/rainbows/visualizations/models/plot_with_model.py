from ....imports import *
from ..utilities import _add_panel_labels

__all__ = ["plot_with_model", "plot_with_model_and_residuals"]


def plot_with_model(
    self,
    quantity="flux",
    errorbar=True,
    text=True,
    ax=None,
    scaling=1,
    data_plotkw={},
    data_errorbarkw={},
    model_plotkw={},
    minimum_acceptable_ok=1e-10,
    filename=None,
    **kw,
):
    """
    Produce a plot of multiple light curves stacked with an offset,
    with models overplotted (like `plot_lightcurves` with model).

    Parameters
    ----------
    quantity : string
        Which quantity to plot? Options include ['flux', 'residuals'].
    errorbar : boolean
        Should we plot errorbars?
    text : boolean
        Should we label each wavelength?
    data_plotkw : dict
        A dictionary of keywords passed to `plt.plot` for the data.
    data_errorbarkw : dict
        A dictionary of keywords passed to `plt.errorbar` for the data.
    model_plotkw : dict
        A dictionary of keywords passed to `plt.plot` for the model.
    **kw : dict
        All additional keywords will be passed along
        to `.plot_lightcurves. Please see the docstring
        for that function for the full set of options.
    """

    # limit to flux or residuals
    assert quantity in ["flux", "residuals", "scaled_residuals_plus_one"]

    # plot the data
    plotkw = dict(marker="o", linewidth=0, markeredgecolor="none", zorder=0)
    plotkw.update(**data_plotkw)

    errorbarkw = data_errorbarkw
    if quantity == "residuals":
        kw.update(quantity="residuals_plus_one")
    elif quantity == "flux":
        kw.update(quantity="flux")
    ax = self.plot_lightcurves(
        errorbar=errorbar,
        text=text,
        ax=ax,
        plotkw=plotkw,
        errorbarkw=errorbarkw,
        minimum_acceptable_ok=minimum_acceptable_ok,
        scaling=scaling,
        **kw,
    )

    # plot the model
    plotkw = dict(marker=None, linewidth=1, zorder=1)
    plotkw.update(**model_plotkw)
    if quantity == "residuals":
        kw.update(quantity="ones")
    elif quantity == "flux":
        kw.update(quantity="model")
    self.plot_lightcurves(
        errorbar=False,
        text=False,
        plotkw=plotkw,
        minimum_acceptable_ok=minimum_acceptable_ok,
        ax=ax,
        **kw,
    )
    if filename is not None:
        self.savefig(filename)
    return ax


def plot_with_model_and_residuals(
    self,
    figsize=(8, 6),
    filename=None,
    histogram=True,
    histogramkw={},
    residual_scaling=1,
    ax=None,
    label_scatter="{measured*1e6:.0f} ppm [{measured/expected:.1f}x]",
    **kw,
):
    """
    Produce a plot of multiple light curves stacked with an offset,
    with models overplotted (like `plot_lightcurves` with model),
    along with residuals in an adjacent panel.

    Parameters
    ----------
    figsize : tuple
        The figure size for the two-panel side-by-side plot.
    histogram : bool
        Should histograms be plotted for the residuals?
    histogramkw : dict
        Keywords to be passed to plot histograms.
    **kw : dict
        All additional keywords will be passed along
        to `.plot_with_models. Please see the docstring
        for that function for the full set of options.
    """
    if histogram:
        N_columns = 3
        gridspec_kw = dict(width_ratios=[1, 1, 0.1])
    else:
        N_columns = 2
        gridspec_kw = dict()

    if ax is None:
        fi, ax = plt.subplots(
            1,
            N_columns,
            figsize=figsize,
            dpi=300,
            sharey=True,
            # sharex=True,
            constrained_layout=True,
            gridspec_kw=gridspec_kw,
        )

    self.plot_with_model(ax=ax[0], **kw)
    spacing = ax[0]._most_recent_chromatic_plot_spacing
    if "spacing" not in kw:
        kw.update(spacing=spacing)
    if label_scatter is not False:
        kw.update(label_scatter=label_scatter)

    self.plot_with_model(
        quantity="residuals",
        scaling=residual_scaling,
        ax=ax[1],
        **kw,
    )
    if residual_scaling != 1:
        ax[1].set_ylabel(
            r"Residuals $\times$ " + f"{residual_scaling}" + " (+ offsets)"
        )
    else:
        ax[1].set_ylabel("Residuals (+ offsets)")

    if histogram:
        kw = dict(
            orientation="horizontal",
            quantity="residuals",
            expected=True,
            color="auto",
            scaling=residual_scaling,
        )
        kw.update(**histogramkw)
        for i in range(self.nwave):
            self.plot_histogram(i, offset=1 - i * spacing, ax=ax[2], **kw)
        plt.axis("off")

    if filename is not None:
        self.savefig(filename)

    return ax
