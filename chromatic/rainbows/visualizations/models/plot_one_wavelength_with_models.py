from ....imports import *
from ..utilities import *

__all__ = ["plot_one_wavelength_with_models", "animate_with_models"]


def _create_title(s):
    """
    A helper function to make better titles.
    """
    if s == "model":
        return "Relative Flux"
    elif "model" in s:
        return "Flux (" + s.replace("_", " ").title().replace(" Model", "") + ")"
    else:
        return s.title()


def plot_one_wavelength_with_models(
    self,
    i_wavelength,
    models=["model", "systematics_model", "planet_model"],
    vlimits_data=[None, None],
    vspan_residuals=None,
    panelsize=(3.25, 2),
    orientation="vertical",
    errorbar=True,
    text=True,
    data_errorbarkw={},
    model_plotkw={},
    wavelength_textkw={},
    cmap=None,
    vmin=None,
    vmax=None,
    w_unit="micron",
    t_unit="day",
    label="inside",
    label_textkw={"color": "black"},
    animation=False,
    minimum_acceptable_ok=1e-10,
    **kw,
):
    """
    Produce a multi-panel plot figure comparing
    the data to model components. The first panel
    will be the data ('flux'), the last panel will
    be the residuals ('flux' - 'model'), and the middle
    panel(s) will be the model(s) connecting them.

    Parameters
    ----------
    i_wavelength : int
        The index of the wavelength to plot.
    models : list, optional
        The list of model quantities to display,
        each in its own panel. A popular option might be:
        `models = ['systematics_model', 'planet_model']
        Whatever you include here does not affect the
        values of the Residuals panel, which are always
        exactly 'flux' - 'model'.
    vlimits_data : list, optional
        The value limits for the data and model colormap,
        specified as `vlimits_data = [vmin, vmax]`.
        If `None` or `[None, None]`, the limits will be
        replaced with the [1,99] percentiles of the flux.
    vspan_residuals : float, optional
        The distance away from zero where the colormap
        for the residuals will cut off. If `None`, the
        span will be set to 3 standard deviations.
    panelsize : tuple, optional
        The size of each panel of the multipanel plot.
        It should probably be wider than it is tall.
    orientation : str, optional
        The direction in which the panels should spread out.
        Options include ['vertical', 'horizontal'].
    errorbar : boolean, optional
        Should we plot errorbars?
    text : boolean, optional
        Should we label each wavelength?
    data_errorbarkw : dict, optional
        A dictionary of keywords passed to `plt.errorbar`
        so you can have more detailed control over the data
        appearance. Common keyword arguments might include:
        `[alpha, elinewidth, color, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
    model_plotkw : dict, optional
        A dictionary of keywords passed to `plt.plot`
        so you can have more detailed control over the model
        appearance. Common keyword arguments might include:
        `[alpha, clip_on, zorder, marker, markersize,
          linewidth, linestyle, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
    wavelength_textkw : dict, optional
        A dictionary passed to `plt.text` for the wavelength label
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    cmap : str, Colormap, optional
        The color map to use for expressing wavelength
    vmin : Quantity, optional
        The minimum value to use for the wavelength colormap
    vmax : Quantity, optional
        The maximum value to use for the wavelength colormap
    w_unit : str, Unit, optional
        The unit for plotting wavelengths.
    t_unit : str, Unit, optional
        The unit for plotting times.
    label : bool, optional
        Should we add (a), (b), (c) labels to the panels,
        and where? False or None give no labels, 'inside'
        sets them in the upper left inside corner, 'outside'
        sets them above the upper left corner.
    label_textkw : dict, optional
        A dictionary passed to `plt.text` for the (a), (b), (c) label
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    animation : bool, optional
        Are we in the middle of an animation
        (so shouldn't make a new figure)?
    **kw : dict, optional
        All additional keywords will be ignored.
    """

    # make sure the units work
    w_unit, t_unit = u.Unit(w_unit), u.Unit(t_unit)

    # make sure we use only models that exist
    models_that_exist = []
    for m in models:
        if self.get(m) is None:
            cheerfully_suggest(f"'{m}' doesn't exist and will be skipped.")
        else:
            models_that_exist.append(m)
    models = models_that_exist

    # make sure vertical limits are set for data and model
    percentiles = [0, 100]
    if vlimits_data == None or vlimits_data == [None, None]:
        vlimits_data = np.percentile(self.flux, percentiles)
    elif vlimits_data[0] == None:
        vlimits_data[0] = np.percentile(self.flux, percentiles[0])
    elif vlimits_data[1] == None:
        vlimits_data[1] = np.percentile(self.flux, percentiles[1])

    # make sure vertical limits are set for the residuals
    if vspan_residuals is None:
        vspan_residuals = np.nanstd(self.residuals) * 5

    # make sure that the wavelength-based colormap is defined
    self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

    # set up the panels
    N = len(models) + 1
    if orientation == "vertical":
        rows, columns = N, 1
        sharex, sharey = "all", "none"
        vspan_residual_ratio = np.maximum(
            vspan_residuals / (vlimits_data[1] - vlimits_data[0]), 0.5
        )
        height_ratios = [1] * (N - 1) + [vspan_residual_ratio]
        gridspec_kw = dict(height_ratios=height_ratios)
        figsize = (panelsize[0], np.sum(height_ratios) * panelsize[1])
    elif orientation == "horizontal":
        rows, columns = 1, N
        sharex, sharey = "all", "all"
        gridspec_kw = dict()
        figsize = (panelsize[0] * N, panelsize[1])
    else:
        raise ValueError(
            f"Please choose an orientation from ['vertical', 'horizontal']."
        )

    # allow animation by making an option to not remake figure and axes
    setup_keys = ["fi", "ax", "vlimits_data", "vspan_residuals"]
    try:
        # try to pull the pre-existing figures and axes (if making animation)
        fi = self._animate_with_models_setup["fi"]
        ax = self._animate_with_models_setup["ax"]
        vlimits_data = self._animate_with_models_setup["vlimits_data"]
        vspan_residuals = self._animate_with_models_setup["vspan_residuals"]
        assert animation == True
        for a in ax:
            a.cla()
    except (AttributeError, AssertionError):
        # otherwise, create  new figure and axes
        fi, ax = plt.subplots(
            rows,
            columns,
            figsize=figsize,
            dpi=300,
            sharex=sharex,
            sharey=sharey,
            gridspec_kw=gridspec_kw,
            constrained_layout=True,
            facecolor="white",
        )
        self._animate_with_models_setup = {}
        for k in setup_keys:
            self._animate_with_models_setup[k] = vars()[k]

    # get the automatic color for this particular wavelength
    w = self.wavelength[i_wavelength]
    color = self.get_wavelength_color(w)

    # set the default xlabel
    xlabel = f"Time ({t_unit.to_string('latex_inline')})"

    # set defaults for the model line plots
    plotkw = dict(marker=None, linewidth=1, zorder=1, color=color, alpha=0.5)
    plotkw.update(**model_plotkw)

    # set default for the data errorbar plots
    default_markersize = plt.matplotlib.rcParams["lines.markersize"]
    errorbarkw = dict(
        color=color,
        linewidth=0,
        elinewidth=int(errorbar),
        marker="o",
        # markersize=int(errorbar == False) * default_markersize,
        markeredgecolor="none",
        alpha=0.5,
        zorder=-1,
    )
    errorbarkw.update(**data_errorbarkw)

    # pull out the raw data and complete model
    data_x, data_y, data_sigma = self.get_ok_data_for_wavelength(
        i_wavelength,
        x="time",
        y="flux",
        sigma="uncertainty",
        minimum_acceptable_ok=minimum_acceptable_ok,
    )
    model_x, model_y, _ = self.get_ok_data_for_wavelength(
        i_wavelength,
        x="time",
        y="model",
        sigma="uncertainty",
        minimum_acceptable_ok=minimum_acceptable_ok,
    )

    # loop through models to include
    for i, (k, a) in enumerate(zip(models + ["residuals"], ax)):

        # choose the model we want to isolate
        if k == "residuals":
            this_model_x, this_model_y, _ = self.get_ok_data_for_wavelength(
                i_wavelength,
                x="time",
                y="ones",
                minimum_acceptable_ok=minimum_acceptable_ok,
            )
        else:
            this_model_x, this_model_y, _ = self.get_ok_data_for_wavelength(
                i_wavelength, x="time", y=k, minimum_acceptable_ok=minimum_acceptable_ok
            )

        # point to the appropriate plotting axes
        plt.sca(a)

        # plot the data
        plt.errorbar(
            data_x.to(t_unit),
            data_y / model_y * this_model_y,
            yerr=data_sigma,
            **errorbarkw,
        )

        # plot the model
        plt.plot(this_model_x.to(t_unit), this_model_y, **plotkw)

        # set the axis labels and limits
        if orientation == "vertical":
            a.set_ylabel(_create_title(k))
            if k == "residuals":
                a.set_xlabel(xlabel)
                a.set_ylim(1 - vspan_residuals, 1 + vspan_residuals)
            else:
                a.set_ylim(*vlimits_data)
        elif orientation == "horizontal":
            a.set_title(_create_title(k))
            if i == 0:
                a.set_ylabel("Relative Flux")
            a.set_ylim(*vlimits_data)
            a.set_xlabel(xlabel)

    # add text labels for the wavelength
    textkw = dict(
        x=0.98,
        y=0.04,
        va="bottom",
        ha="right",
        transform=ax[0].transAxes,
        color=color,
    )
    textkw.update(**wavelength_textkw)
    textkw.update(s=f"{w.to_value(w_unit):.2f} {w_unit.to_string('latex_inline')}")
    if text:
        ax[0].text(**textkw)

    # add (a), (b), (c) labels to the panels, if desired
    if label:
        _add_panel_labels(ax[:], preset=label, **label_textkw)

    # set xlimits
    ax[0].set_xlim(self.time.to_value(t_unit).min(), self.time.to_value(t_unit).max())
    plt.sca(ax[0])


def animate_with_models(
    self,
    filename="animated-lightcurves-with-models.gif",
    fps=5,
    bitrate=None,
    dpi=None,
    orientation="horizontal",
    **kw,
):
    """
    Animate light curves with models by flipping through wavelengths.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif or .html files.
    fps : float
        frames/second of animation
    bitrate : None
        How many bits/second should be allowed in the file?
    dpi : float
        How many dots per inch should go into the animation rendering?
    orientation : str, optional
        The direction in which the panels should spread out.
        Options include ['vertical', 'horizontal'].
    **kw : dict
        All other keywords will be passed to `.plot_one_wavelength_with_models`
    """
    try:
        del self._animate_with_models_setup
    except AttributeError:
        pass

    self.plot_one_wavelength_with_models(
        0, animation=False, orientation=orientation, **kw
    )

    filename = self._label_plot_file(filename)

    # initialize the animator
    writer, displayer = get_animation_writer_and_displayer(
        filename=filename, fps=fps, bitrate=bitrate
    )

    # set up to save frames directly into the animation
    figure = self._animate_with_models_setup["fi"]
    with writer.saving(figure, filename, dpi or figure.get_dpi()):

        # loop over exposures
        for i in tqdm(range(self.nwave), leave=False):

            self.plot_one_wavelength_with_models(
                i, animation=True, orientation=orientation, **kw
            )
            # save this snapshot to a movie frame
            writer.grab_frame()

    # close the figure that was created
    plt.close(figure)

    # display the animation
    from IPython.display import display

    try:
        display(displayer(filename, embed=True))
    except TypeError:
        display(displayer(filename))
