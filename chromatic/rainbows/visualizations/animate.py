from ...imports import *
from .utilities import *

__all__ = [
    "_setup_animated_scatter",
    "_setup_animate_spectra",
    "animate_spectra",
    "_setup_animate_lightcurves",
    "animate_lightcurves",
]


def _setup_animated_scatter(self, ax=None, figurekw={}, scatterkw={}, textkw={}):
    """
    Wrapper to set up the basics of animate-able plot.

    This works for any general plot that has a single-color
    line or set of points (using `plt.plot`), with a text
    label in the upper right corner.

    Parameters
    ----------
    ax : Axes, optional
        The axes into which the plot should be drawn.
        If None, a new one will be created.
    figurekw : dict, optional
        A dictionary of keywords to be passed to `plt.figure`
    scatterkw : dict, optional
        A dictionary of keywords to be passed to `plt.scatter`
    textkw : dict, optional
        A dictionary of keywords to be passed to `plt.text`
    """

    # make sure the ax and figure are defined
    if ax is None:
        kw = dict(facecolor="white")
        kw.update(**figurekw)
        fig, ax = plt.subplots(**kw)
    else:
        fig = ax.get_figure()

    plt.sca(ax)

    this_scatterkw = dict(c=[], cmap=self.cmap, norm=self.norm)
    this_scatterkw.update(**scatterkw)
    scatter = plt.scatter([], [], **this_scatterkw)

    this_textkw = dict(
        x=0.98, y=0.96, s="", ha="right", va="top", transform=ax.transAxes
    )
    this_textkw.update(**textkw)
    text = plt.text(**this_textkw)

    plt.title(self.get("title"))

    # return a dictionary with things that will be useful to hang onto
    return dict(fi=fig, ax=ax, scatter=scatter, text=text)


def _setup_animate_lightcurves(
    self,
    ax=None,
    quantity="flux",
    xlim=[None, None],
    ylim=[None, None],
    cmap=None,
    vmin=None,
    vmax=None,
    ylabel=None,
    scatterkw={},
    textkw={},
):
    """
    Setup an animation to how the lightcurve changes
    as we flip through every wavelength.

    Parameters
    ----------
    ax : Axes, optional
        The axes into which this animated plot should go.
    quantity : string, optional
        Which fluxlike quantity should be retrieved? (default = 'flux')
    xlim : tuple, optional
        Custom xlimits for the plot
    ylim : tuple, optional
        Custom ylimits for the plot
    cmap : str, Colormap, optional
        The color map to use for expressing wavelength
    vmin : Quantity, optional
        The minimum value to use for the wavelength colormap
    vmax : Quantity, optional
        The maximum value to use for the wavelength colormap
    scatterkw : dict, optional
        A dictionary of keywords to be passed to `plt.scatter`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[s, c, marker, alpha, linewidths, edgecolors, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
    textkw : dict
        A dictionary of keywords passed to `plt.text`
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    if (quantity == "flux") and (self._is_probably_normalized() == False):
        cheerfully_suggest(
            f"""
        It's not 100% obvious that {self} has been normalized.
        If you're expecting your animation to wobble near 1,
        please try normalizing before animating, perhaps with
        the `.normalize()` action. If not, please ignore this!
        """
        )

    self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

    with quantity_support():

        # keep track of the things needed for the animation
        self._animate_lightcurves_components = self._setup_animated_scatter(
            ax=ax, scatterkw=scatterkw, textkw=textkw
        )

        ax = self._animate_lightcurves_components["ax"]

        # set the plot limits
        ax.set_xlim(xlim[0] or np.nanmin(self.time), xlim[1] or np.nanmax(self.time))
        ax.set_ylim(
            ylim[0] or 0.995 * np.nanmin(self.get(quantity)),
            ylim[1] or 1.005 * np.nanmax(self.get(quantity)),
        )
        # set the axis labels
        ax.set_xlabel(
            f"{self._time_label} ({self.time.unit.to_string('latex_inline')})"
        )
        ax.set_ylabel(ylabel or quantity)

        # guess a good number of digits to round
        ndigits = np.minimum(
            int(np.floor(np.log10(np.min(np.diff(self.wavelength)).value))), 0
        )
        format_code = f"{{:.{-ndigits}f}}"

        def update(frame):
            """
            This function will be called to update each frame
            of the animation.

            Parameters
            ----------
            frame : int
                An integer that will advance with each frame.
            """

            # pull out the x and y values to plot
            x, y, _ = self.get_ok_data_for_wavelength(frame, y=quantity)
            # x = self.time
            # y = self.flux[frame]
            c = self.wavelength[frame].to("micron").value * np.ones(self.ntime)

            # update the label in the corner
            self._animate_lightcurves_components["text"].set_text(
                f"w = {format_code.format(self.wavelength[frame].value)} {self.wavelength.unit.to_string('latex')}"
            )

            # update the plot data
            self._animate_lightcurves_components["scatter"].set_offsets(
                np.transpose([x, y])
            )
            self._animate_lightcurves_components["scatter"].set_array(c)

            return (
                self._animate_lightcurves_components["text"],
                self._animate_lightcurves_components["scatter"],
            )

        # hold onto this update function in case we need it elsewhere
        self._animate_lightcurves_components["update"] = update


def animate_lightcurves(
    self,
    filename="animated-lightcurves.gif",
    fps=5,
    dpi=None,
    bitrate=None,
    **kwargs,
):
    """
    Create an animation to show how the lightcurve changes
    as we flip through every wavelength.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif or .html files.
    fps : float
        frames/second of animation
    ax : Axes
        The axes into which this animated plot should go.
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    cmap : str,
        The color map to use for expressing wavelength
    vmin : Quantity
        The minimum value to use for the wavelength colormap
    vmax : Quantity
        The maximum value to use for the wavelength colormap
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[s, c, marker, alpha, linewidths, edgecolors, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
    textkw : dict
        A dictionary of keywords passed to `plt.text`
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    self._setup_animate_lightcurves(**kwargs)

    filename = self._label_plot_file(filename)

    # initialize the animator
    writer, displayer = get_animation_writer_and_displayer(
        filename=filename, fps=fps, bitrate=bitrate
    )

    # set up to save frames directly into the animation
    figure = self._animate_lightcurves_components["fi"]
    with writer.saving(figure, filename, dpi or figure.get_dpi()):
        for i in tqdm(range(self.nwave), leave=False):
            self._animate_lightcurves_components["update"](i)
            writer.grab_frame()

    # close the figure that was created
    plt.close(figure)

    # display the animation
    from IPython.display import display

    try:
        display(displayer(filename, embed=True))
    except TypeError:
        display(displayer(filename))


def _setup_animate_spectra(
    self,
    ax=None,
    quantity="flux",
    xlim=[None, None],
    ylim=[None, None],
    cmap=None,
    vmin=None,
    vmax=None,
    ylabel=None,
    scatterkw={},
    textkw={},
):
    """
    Setup an animation to show how the spectrum changes
    as we flip through every timepoint.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif files.
    ax : Axes
        The axes into which this animated plot should go.
    quantity : string
        Which fluxlike quantity should be retrieved? (default = 'flux')
    fps : float
        frames/second of animation
    figsize : tuple
        (width, height) of the figure
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    cmap : str,
        The color map to use for expressing wavelength
    vmin : Quantity
        The minimum value to use for the wavelength colormap
    vmax : Quantity
        The maximum value to use for the wavelength colormap
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[s, c, marker, alpha, linewidths, edgecolors, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
    textkw : dict
        A dictionary of keywords passed to `plt.text`
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """

    if (quantity == "flux") and (self._is_probably_normalized() == False):
        cheerfully_suggest(
            f"""
        It's not 100% obvious that {self} has been normalized.
        If you're expecting your animation to wobble near 1,
        please try normalizing before animating, perhaps with
        the `.normalize()` action. If not, please ignore this!
        """
        )
    self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

    with quantity_support():

        # keep track of the things needed for the animation
        self._animate_spectra_components = self._setup_animated_scatter(
            ax=ax, scatterkw=scatterkw, textkw=textkw
        )

        ax = self._animate_spectra_components["ax"]

        """if self.wscale == "log":
            plt.xscale("log")
            formatter = plt.matplotlib.ticker.StrMethodFormatter("{x:.1g}")
            ax.xaxis.set_major_formatter(formatter)
            ax.xaxis.set_minor_formatter(formatter)"""

        # set the plot limits
        ax.set_xlim(
            xlim[0] or np.nanmin(self.wavelength), xlim[1] or np.nanmax(self.wavelength)
        )
        ax.set_ylim(
            ylim[0] or 0.995 * np.nanmin(self.get(quantity)),
            ylim[1] or 1.005 * np.nanmax(self.get(quantity)),
        )
        # set the axis labels
        ax.set_xlabel(
            f"{self._wave_label}  ({self.wavelength.unit.to_string('latex_inline')})"
        )
        ax.set_ylabel(ylabel or quantity)

        # guess a good number of digits to round
        ndigits = np.minimum(
            int(np.floor(np.log10(np.min(np.diff(self.time)).value))), 0
        )
        format_code = f"{{:.{-ndigits}f}}"

        def update(frame):
            """
            This function will be called to update each frame
            of the animation.

            Parameters
            ----------
            frame : int
                An integer that will advance with each frame.
            """

            # pull out the x and y values to plot
            x, y, _ = self.get_ok_data_for_time(frame, y=quantity)
            # x = self.wavelength
            # y = self.flux[:, frame]
            c = self.wavelength.to("micron").value

            # update the label in the corner
            self._animate_spectra_components["text"].set_text(
                f"t = {format_code.format(self.time[frame].value)} {self.time.unit.to_string('latex')}"
            )

            # update the plot data
            self._animate_spectra_components["scatter"].set_offsets(
                np.transpose([x, y])
            )
            self._animate_spectra_components["scatter"].set_array(c)

            return (
                self._animate_spectra_components["text"],
                self._animate_spectra_components["scatter"],
            )

        # hold onto this update function in case we need it elsewhere
        self._animate_spectra_components["update"] = update


def animate_spectra(
    self, filename="animated-spectra.gif", fps=5, dpi=None, bitrate=None, **kw
):
    """
    Create an animation to show how the spectrum changes
    as we flip through every timepoint.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif files.
    ax : Axes
        The axes into which this animated plot should go.
    fps : float
        frames/second of animation
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    cmap : str,
        The color map to use for expressing wavelength
    vmin : Quantity
        The minimum value to use for the wavelength colormap
    vmax : Quantity
        The maximum value to use for the wavelength colormap
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
        so you can have more detailed control over the plot
        appearance. Common keyword arguments might include:
        `[s, c, marker, alpha, linewidths, edgecolors, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
    textkw : dict
        A dictionary of keywords passed to `plt.text`
        so you can have more detailed control over the text
        appearance. Common keyword arguments might include:
        `[alpha, backgroundcolor, color, fontfamily, fontsize,
          fontstyle, fontweight, rotation, zorder]` (and more)
        More details are available at
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """

    self._setup_animate_spectra(**kw)

    filename = self._label_plot_file(filename)
    # initialize the animator
    writer, displayer = get_animation_writer_and_displayer(
        filename=filename, fps=fps, bitrate=bitrate
    )

    # set up to save frames directly into the animation
    figure = self._animate_spectra_components["fi"]
    with writer.saving(figure, filename, dpi or figure.get_dpi()):
        for i in tqdm(range(self.ntime), leave=False):
            self._animate_spectra_components["update"](i)
            writer.grab_frame()

    # close the figure that was created
    plt.close(figure)

    # display the animation
    from IPython.display import display

    display(displayer(filename))
