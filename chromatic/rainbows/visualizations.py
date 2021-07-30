from ..imports import *


def imshow(
    self,
    ax=None,
    quantity="flux",
    w_unit="micron",
    t_unit="hour",
    aspect="auto",
    colorbar=True,
    origin="upper",
    **kw,
):
    """
    imshow flux as a function of time (x = time, y = wavelength, color = flux).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    w_unit : str, astropy.unit.Unit
        The unit for plotting wavelengths.
    t_unit : str, astropy.unit.Unit
        The unit for plotting times.

    """

    # self.speak(f'imshowing')
    if ax is None:
        ax = plt.subplot()

    w_unit, t_unit = u.Unit(w_unit), u.Unit(t_unit)

    if self.wscale == "linear":
        extent = [
            (min(self.time) / t_unit).decompose(),
            (max(self.time) / t_unit).decompose(),
            (max(self.wavelength) / w_unit).decompose(),
            (min(self.wavelength) / w_unit).decompose(),
        ]
    elif self.wscale == "log":
        extent = [
            (min(self.time) / t_unit).decompose(),
            (max(self.time) / t_unit).decompose(),
            np.log10(max(self.wavelength) / w_unit),
            np.log10(min(self.wavelength) / w_unit),
        ]
    else:
        raise RuntimeError("Can't imshow without knowing wscale.")

    with quantity_support():
        plt.sca(ax)
        plt.imshow(
            self.fluxlike[quantity],
            extent=extent,
            aspect=aspect,
            origin=origin,
            **kw,
        )
        if self.wscale == "linear":
            plt.ylabel(f"Wavelength ({w_unit.to_string('latex_inline')})")
        elif self.wscale == "log":
            plt.ylabel(
                r"log$_{10}$" + f"[Wavelength/({w_unit.to_string('latex_inline')})]"
            )
        plt.xlabel(f"Time ({t_unit.to_string('latex_inline')})")
        if colorbar:
            plt.colorbar(ax=ax)
    return ax


def plot(self, ax=None, spacing=None, plotkw={}, fontkw={}):
    """
    Plot flux as sequence of offset light curves.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    spacing : None, float
        The spacing between light curves.
        (Might still change how this works.)
    plowkw : dict
        A dictionary of keywords passed to `plt.plot`
    textkw : dict
        A dictionary of keywords passed to `plt.text`
    """
    from chromatic import viz

    viz.wavelength_plot(
        self.flux,
        self.time,
        self.wavelength,
        step_size=spacing,
        ax=ax,
        plotkw=plowkw,
        fontkw=fontkw,
    )


def _setup_animated_plot(self, ax=None, figsize=None, plotkw={}, textkw={}):
    """
    Wrapper to set up the basics of animate-able plot.

    This works for any general plot that has a single-color
    line or set of points (using `plt.plot`), with a text
    label in the upper right corner.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes into which the plot should be drawn.
        If None, a new one will be created.
    figsize : tuple
        (width, height) of the figure, if one needs
        to be created (= if ax isn't specified)
    plotkw : dict
        A dictionary of keywords to be passed to `plt.plot`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    """

    # make sure the ax and figure are defined
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    default_plotkw = dict()
    plot = plt.plot([], [], **default_plotkw, **plotkw)[0]

    default_textkw = dict(
        x=0.98, y=0.96, s="", ha="right", va="top", transform=ax.transAxes
    )
    text = plt.text(**default_textkw, **textkw)

    # return a dictionary with things that will be useful to hang onto
    return dict(fig=fig, ax=ax, plot=plot, text=text)


def animate_lightcurves(
    self,
    filename="animated-lightcurves.gif",
    fps=5,
    figsize=None,
    xlim=[None, None],
    ylim=[None, None],
    plotkw={},
    textkw={},
):
    """
    Create an animation that shows how the lightcurve changes
    as we flip through every wavelength.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif files.
    fps : float
        frames/second of animation
    figsize : tuple
        (width, height) of the figure
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    plotkw : dict
        A dictionary of keywords to be passed to `plt.plot`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    """

    with quantity_support():

        # keep track of the things needed for the animation
        self._animate_lightcurves_components = self._setup_animated_plot(
            figsize=figsize, plotkw=plotkw, textkw=textkw
        )

        ax = self._animate_lightcurves_components["ax"]

        # set the plot limits
        ax.set_xlim(xlim[0] or np.min(self.time), xlim[1] or np.max(self.time))
        ax.set_ylim(
            ylim[0] or 0.995 * np.min(self.flux),
            ylim[1] or 1.005 * np.max(self.flux),
        )
        # set the axis labels
        ax.set_xlabel(f"Time ({self.time.unit.to_string('latex_inline')})")
        ax.set_ylabel(f"Relative Flux")

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
            x = self.time
            y = self.flux[frame]

            # update the label in the corner
            self._animate_lightcurves_components["text"].set_text(
                f"w = {self.wavelength[frame].value:0.2f} {self.wavelength.unit.to_string('latex')}"
            )

            # update the plot data
            self._animate_lightcurves_components["plot"].set_data(x, y)

            return (
                self._animate_lightcurves_components["text"],
                self._animate_lightcurves_components["plot"],
            )

        # hold onto this update function in case we need it elsewhere
        self._animate_lightcurves_components["update"] = update

        # make and save the animation
        ani = matplotlib.animation.FuncAnimation(
            self._animate_lightcurves_components["fig"],
            update,
            frames=np.arange(0, self.nwave),
            blit=True,
        )
        ani.save(filename, fps=fps)


def animate_spectra(
    self,
    filename="animated-spectra.gif",
    fps=5,
    figsize=None,
    xlim=[None, None],
    ylim=[None, None],
    plotkw={},
    textkw={},
):
    """
    Create an animation that shows how the lightcurve changes
    as we flip through every wavelength.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif files.
    fps : float
        frames/second of animation
    figsize : tuple
        (width, height) of the figure
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    plotkw : dict
        A dictionary of keywords to be passed to `plt.plot`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    """

    with quantity_support():

        # keep track of the things needed for the animation
        self._animate_spectra_components = self._setup_animated_plot(
            figsize=figsize, plotkw=plotkw, textkw=textkw
        )

        ax = self._animate_spectra_components["ax"]

        # set the plot limits
        ax.set_xlim(
            xlim[0] or np.min(self.wavelength), xlim[1] or np.max(self.wavelength)
        )
        ax.set_ylim(
            ylim[0] or 0.995 * np.min(self.flux),
            ylim[1] or 1.005 * np.max(self.flux),
        )
        # set the axis labels
        ax.set_xlabel(f"Wavelength ({self.wavelength.unit.to_string('latex_inline')})")
        ax.set_ylabel(f"Relative Flux")

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
            x = self.wavelength
            y = self.flux[:, frame]

            # update the label in the corner
            self._animate_spectra_components["text"].set_text(
                f"t = {self.time[frame].value:0.2f} {self.time.unit.to_string('latex')}"
            )

            # update the plot data
            self._animate_spectra_components["plot"].set_data(x, y)

            return (
                self._animate_spectra_components["text"],
                self._animate_spectra_components["plot"],
            )

        # hold onto this update function in case we need it elsewhere
        self._animate_spectra_components["update"] = update

        # make and save the animation
        ani = matplotlib.animation.FuncAnimation(
            self._animate_spectra_components["fig"],
            update,
            frames=np.arange(0, self.ntime),
            blit=True,
        )
        ani.save(filename, fps=fps)
