from ..imports import *


def _setup_wavelength_colors(self, cmap=None, vmin=None, vmax=None):
    """
    Set up a color map and normalization function for
    colors datapoints by their wavelengths.

    Parameters
    ----------
    cmap : str, matplotlib.colors.Colormap
        The color map to use.
    vmin : astropy.units.Quantity
        The wavelength at the bottom of the cmap.
    vmax : astropy.units.Quantity
        The wavelength at the top of the cmap.
    """

    # populate the cmap object
    self.cmap = cm.get_cmap(cmap)

    vmin = vmin or np.nanmin(self.wavelength)
    vmax = vmax or np.nanmax(self.wavelength)

    if self.wscale in ["?", "linear"]:
        self.norm = col.Normalize(
            vmin=vmin.to("micron").value, vmax=vmax.to("micron").value
        )
    elif self.wscale in ["log"]:
        self.norm = col.LogNorm(
            vmin=vmin.to("micron").value, vmax=vmax.to("micron").value
        )


def _make_sure_cmap_is_defined(self, cmap=None, vmin=None, vmax=None):
    """
    A helper function that can be called at the start of
    any plot that that's using wavelength-colors to make
    sure that the wavelength-based colormap has been
    defined.
    """

    if hasattr(self, "cmap"):
        if (cmap is not None) or (vmin is not None) or (vmax is not None):
            warnings.warn(
                """
            It looks like you're trying to set up a new custom
            cmap and/or wavelength normalization scheme. You
            should be aware that a cmap has already been defined
            for this object; if you're visualizing the same
            rainbow in different ways, we strongly suggest
            that you not change the cmap or normalization
            between them, for visual consistency.
            """
            )
        else:
            return
    self._setup_wavelength_colors(cmap=cmap, vmin=vmin, vmax=vmax)


def get_wavelength_color(self, wavelength):
    """
    Determine the color corresponding to one or more wavelengths.

    Parameters
    ----------
    wavelength : astropy.units.Quantity
        The wavelength value(s), either an individual
        wavelength or an array of N wavelengths.

    Returns
    -------
    colors : np.array
        An array of RGBA colors [or an (N,4) array].
    """
    w_unitless = wavelength.to("micron").value
    normalized_w = self.norm(w_unitless)
    return self.cmap(normalized_w)


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


def plot(
    self, ax=None, spacing=None, cmap=None, vmin=None, vmax=None, plotkw={}, textkw={}
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
    plowkw : dict
        A dictionary of keywords passed to `plt.plot`
    textkw : dict
        A dictionary of keywords passed to `plt.text`
    """

    # make sure that the wavelength-based colormap is defined
    self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

    time_unit = self.time.unit.to_string("latex_inline")
    wave_unit = self.wavelength.unit.to_string("latex_inline")

    min_time = np.min(self.time)

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
    warnings.warn(
        """
    It's not clear if/how this object has been normalized.
    Be aware that the baseline flux levels may therefore
    be a little bit funny in .plot()."""
    )
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
                plt.plot(self.time, plot_flux, **this_plotkw)

                # add text labels next to each light curve
                this_textkw = dict(va="center", color=color)
                this_textkw.update(**textkw)
                plt.annotate(
                    f"{w.value:4.2f} {wave_unit}",
                    (min_time, 1 - (i + 0.5) * spacing),
                    **this_textkw,
                )

        # add text labels to the plot
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
