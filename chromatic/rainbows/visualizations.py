from ..imports import *


def name2color(name):
    """
    Return the 3-element RGB array of a given color name.

    Parameters
    ----------
    name : str
        The name of a color

    Returns
    -------
    rgb : tuple
        3-element RGB color
    """
    try:
        color_hex = col.cnames[name]
        return col.hex2color(color_hex)
    except KeyError:
        warnings.warn(f"The color {name} can't be found. (Returning black.)")
        return (0.0, 0.0, 0.0)


def one2another(bottom="white", top="red", alphabottom=1.0, alphatop=1.0, N=256):
    """
    Create a cmap that goes smoothly (linearly in RGBA) from "bottom" to "top".
    """
    rgb_bottom, rgb_top = name2color(bottom), name2color(top)
    r = np.linspace(rgb_bottom[0], rgb_top[0], N)
    g = np.linspace(rgb_bottom[1], rgb_top[1], N)
    b = np.linspace(rgb_bottom[2], rgb_top[2], N)
    a = np.linspace(alphabottom, alphatop, N)
    colors = np.transpose(np.vstack([r, g, b, a]))
    cmap = co.ListedColormap(colors, name="{bottom}2{top}".format(**locals()))
    return cmap


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
        None uses half the standard dev of entire flux data.
    plowkw : dict
        A dictionary of keywords passed to `plt.plot`
    textkw : dict
        A dictionary of keywords passed to `plt.text`
    """

    time = self.time.value
    time_unit = self.time.unit.to_string("latex_inline")
    wavelength = self.wavelength.value
    waveunit = self.wavelength.unit.to_string("latex_inline")

    fcadences = 0.05
    start_pcad = int(np.ceil(len(time) * fcadences))
    end_pcad = int(np.floor(len(time) * fcadences))

    cont_time = np.zeros(len(time))
    cont_time[0:start_pcad] = 1
    cont_time[end_pcad:-1] = 1

    min_time = np.min(time)

    if spacing is None:
        spacing = 3 * np.nanstd(self.flux)

    nsteps = len(wavelength)

    if ax is None:
        ax = plt.subplot()
    plt.sca(ax)
    # plt.figure(figsize = (8,12) )

    for i, w in enumerate(wavelength):

        lc = self.flux[i, :]

        if np.any(np.isfinite(lc)):
            # normalize this light curve to one.
            cont_level = np.nanmedian(lc[cont_time == 1])
            plot_flux = -i * spacing + lc  # / cont_level

            color = (0, 0.3 - 0.3 * (i / nsteps), i / nsteps)
            assert (np.isfinite(plot_flux) == True).all()

            plt.plot(time, plot_flux, ".--", color=color, **plotkw)
            plt.annotate(
                f"{wavelength[i]:4.2f} {waveunit}",
                (min_time, 1 - i * spacing - 0.5 * spacing),
                color=color,
                **fontkw,
            )

    plt.xlabel(f"Time ({time_unit})")
    plt.ylabel("Relative Flux (+ offsets)")


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
        animator = ani.FuncAnimation(
            self._animate_lightcurves_components["fig"],
            update,
            frames=np.arange(0, self.nwave),
            blit=True,
        )
        animator.save(filename, fps=fps)


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
        animator = ani.FuncAnimation(
            self._animate_spectra_components["fig"],
            update,
            frames=np.arange(0, self.ntime),
            blit=True,
        )
        animator.save(filename, fps=fps)
