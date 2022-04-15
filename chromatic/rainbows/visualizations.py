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

    Parameters
    ----------
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength.
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap.
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap.
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
    t_unit="day",
    colorbar=True,
    aspect="auto",
    **kw,
):
    """
    imshow flux as a function of time (x = time, y = wavelength, color = flux).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    quantity : str
        The fluxlike quantity to imshow.
        (Must be a key of `rainbow.fluxlike`).
    w_unit : str, astropy.unit.Unit
        The unit for plotting wavelengths.
    t_unit : str, astropy.unit.Unit
        The unit for plotting times.
    colorbar : bool
        Should we include a colorbar?
    aspect : str
        What aspect ratio should be used for the imshow?
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
        ylabel = f"Wavelength ({w_unit.to_string('latex_inline')})"

    elif self.wscale == "log":
        extent = [
            (min(self.time) / t_unit).decompose(),
            (max(self.time) / t_unit).decompose(),
            np.log10(max(self.wavelength) / w_unit),
            np.log10(min(self.wavelength) / w_unit),
        ]
        ylabel = r"log$_{10}$" + f"[Wavelength/({w_unit.to_string('latex_inline')})]"

    else:
        message = f"""
        The wavelength scale for this rainbow is '{self.wscale}'.
        It's hard to imshow something with a wavelength axis
        that isn't uniform in linear or logarithmic space, so
        we're giving up and just using the wavelength index
        as the wavelength axis.

        If you want a real wavelength axis, one solution would
        be to bin your wavelengths to a more uniform grid with
        `rainbow.bin(R=...)` (for logarithmic wavelengths) or
        `rainbow.bin(dw=...)` (for linear wavelengths)
        """
        warnings.warn(message)
        extent = [
            (min(self.time) / t_unit).decompose(),
            (max(self.time) / t_unit).decompose(),
            self.nwave,
            0,
        ]
        ylabel = "Wavelength Index"

    # define some default keywords
    imshow_kw = dict(interpolation="nearest")
    imshow_kw.update(**kw)
    with quantity_support():
        plt.sca(ax)
        plt.imshow(
            self.fluxlike[quantity],
            extent=extent,
            aspect=aspect,
            origin="upper",
            **imshow_kw,
        )
        plt.ylabel(ylabel)
        plt.xlabel(f"Time ({t_unit.to_string('latex_inline')})")
        if colorbar:
            plt.colorbar(
                ax=ax,
                label=u.Quantity(self.fluxlike[quantity]).unit.to_string(
                    "latex_inline"
                ),
            )
    return ax


def imshow_fluxlike_quantities(
    self, quantities=None, maxcol=3, panel_size=(5, 4), **kw
):
    """
    imshow fluxlikes as a function of time (x = time, y = wavelength, color = flux).

    Parameters
    ----------
    quantities : list like
        The fluxlike quantity to imshow.
    maxcol : int
        The maximum number of columns to show (Optional).
    panel_size : tuple
        The (approximate) size of a single panel, which will
        be used to set the overall figsize based on the
        number of rows and columns (Optional).
    kw : dict
        Additional keywords will be passed on to `imshow`

    """

    # decide which quantities to plot
    if quantities is None:
        allkeys = self.fluxlike.keys()
    else:
        allkeys = quantities[:]

    # set up the geometry of the grid
    if len(allkeys) > maxcol:
        rows = int(np.ceil(len(allkeys) / maxcol))
        cols = maxcol
    else:
        rows = 1
        cols = np.min([len(allkeys), maxcol])

    # create the figure and grid of axes
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(cols * panel_size[0], rows * panel_size[1]),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )

    # make the axes easier to index
    if len(allkeys) > 1:
        ax = axes.flatten()
    else:
        ax = [axes]

    # display each quantity
    for k, key in enumerate(allkeys):
        # make the imshow (or an empty box)
        if key in self.fluxlike.keys():
            self.imshow(quantity=key, ax=ax[k], **kw)
        else:
            ax[k].text(
                0.5,
                0.5,
                f"No {key}",
                transform=ax[k].transAxes,
                ha="center",
                va="center",
            )
        # add a title for each box
        ax[k].set_title(key)

        # hide xlabel except on the bottom row
        if k < (len(allkeys) - cols):
            ax[k].set_xlabel("")
        else:
            ax[k].tick_params(labelbottom=True)

        # hide ylabel except on the left column
        if (k % cols) > 0:
            ax[k].set_ylabel("")

    # hide any additional axes
    if k + 1 <= len(ax):
        for axi in ax[k + 1 :]:
            axi.axis("Off")


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


def plot_quantities(
    self,
    quantities=None,
    data_like="time",
    maxcol=1,
    panel_size=(6, 2),
    x_axis="index",
    **kw,
):
    """
    Plot {data_like}-like quantities as a function of {data_like} index
    or any other {data_like} quantity (such as "time" or "wavelength").

    Parameters
    ----------
    quantities : list like
        The X-like quantity to plot.
    data_like : string
        Whether the quantities are alike to 'time' or 'wave'. Default is 'time'. (Optional)
    maxcol : int
        The maximum number of columns to show (Optional).
    panel_size : tuple
        The size in inches dedicated to each panel, default is (6,2). (Optional)
    x_axis : string
        The quantity to plot on the x_axis, default is index. (Optional)

    """
    # decide which dictionary to plot
    if data_like not in ["time", "wave", "wavelength"]:
        raise Exception("Unknown data_like. Choose from [time, wave]")
    elif data_like == "time":
        like_dict = self.timelike
    else:
        like_dict = self.wavelike

    # decide which quantities to plot
    if quantities is None:
        allkeys = like_dict.keys()
    else:
        allkeys = quantities[:]

    # set up the geometry of the grid
    if len(allkeys) > maxcol:
        rows = int(np.ceil(len(allkeys) / maxcol))
        cols = maxcol
    else:
        rows = 1
        cols = np.min([len(allkeys), maxcol])

    # create the figure and grid of axes
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(cols * panel_size[0], rows * panel_size[1]),
        sharex=True,
        constrained_layout=True,
    )
    # make the axes easier to index
    if len(allkeys) > 1:
        ax = axes.flatten()
    else:
        ax = [axes]

    # set x_axis variable
    if x_axis.lower() == "index":
        xaxis = np.arange(0, len(like_dict[list(like_dict.keys())[0]]))
        xlab = f"{data_like} Index"
    else:
        if x_axis in like_dict.keys():
            xaxis = like_dict[x_axis]
            xlab = x_axis
        else:
            raise Exception("Desired x_axis quantity is not in given dictionary")

    # display each quantity
    for k, key in enumerate(allkeys):
        # make the plot (or an empty box)
        if key in like_dict.keys():
            ax[k].plot(
                xaxis, like_dict[key], color=plt.cm.viridis(k / len(allkeys)), **kw
            )
            ax[k].set_xlabel(xlab)
            ax[k].set_ylabel(
                f"{key} ({u.Quantity(like_dict[key]).unit.to_string('latex_inline')})"
            )
        else:
            ax[k].text(
                0.5,
                0.5,
                f"No {key}",
                transform=ax[k].transAxes,
                ha="center",
                va="center",
            )

        # add a title for each box
        ax[k].set_title(key)

        # hide xlabel except on the bottom row
        if k < (len(allkeys) - cols):
            ax[k].set_xlabel("")
        else:
            ax[k].tick_params(labelbottom=True)

    # hide any additional axes
    if k + 1 <= len(ax):
        for axi in ax[k + 1 :]:
            axi.axis("Off")


def _setup_animated_scatter(self, ax=None, scatterkw={}, textkw={}):
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
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    """

    # make sure the ax and figure are defined
    if ax is None:
        fig, ax = plt.subplots()
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

    # return a dictionary with things that will be useful to hang onto
    return dict(fig=fig, ax=ax, scatter=scatter, text=text)


def _setup_animate_lightcurves(
    self,
    ax=None,
    xlim=[None, None],
    ylim=[None, None],
    cmap=None,
    vmin=None,
    vmax=None,
    scatterkw={},
    textkw={},
):
    """
    Setup an animation to how the lightcurve changes
    as we flip through every wavelength.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif files.
    fps : float
        frames/second of animation
    ax : matplotlib.axes.Axes
        The axes into which this animated plot should go.
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    """

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
            ylim[0] or 0.995 * np.nanmin(self.flux),
            ylim[1] or 1.005 * np.nanmax(self.flux),
        )
        # set the axis labels
        ax.set_xlabel(f"Time ({self.time.unit.to_string('latex_inline')})")
        ax.set_ylabel(f"Relative Flux")

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
            x = self.time
            y = self.flux[frame]
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
    self, filename="animated-lightcurves.gif", fps=None, dpi=None, **kwargs
):
    """
    Create an animation to show how the lightcurve changes
    as we flip through every wavelength.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif files.
    fps : float
        frames/second of animation
    ax : matplotlib.axes.Axes
        The axes into which this animated plot should go.
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    """
    self._setup_animate_lightcurves(**kwargs)

    # make and save the animation
    animator = ani.FuncAnimation(
        self._animate_lightcurves_components["fig"],
        self._animate_lightcurves_components["update"],
        frames=np.arange(0, self.nwave),
        blit=True,
    )
    animator.save(filename, fps=fps, dpi=dpi, savefig_kwargs=dict(facecolor="white"))
    plt.close()


def _setup_animate_spectra(
    self,
    ax=None,
    xlim=[None, None],
    ylim=[None, None],
    cmap=None,
    vmin=None,
    vmax=None,
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
    ax : matplotlib.axes.Axes
        The axes into which this animated plot should go.
    fps : float
        frames/second of animation
    figsize : tuple
        (width, height) of the figure
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    fps : int
        Frames per second for the animation.
    dpi : float, default: :rc:`savefig.dpi`
        Dots per inch for the movie.
    """

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
            ylim[0] or 0.995 * np.nanmin(self.flux),
            ylim[1] or 1.005 * np.nanmax(self.flux),
        )
        # set the axis labels
        ax.set_xlabel(f"Wavelength ({self.wavelength.unit.to_string('latex_inline')})")
        ax.set_ylabel(f"Relative Flux")

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
            x = self.wavelength
            y = self.flux[:, frame]
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
    self, filename="animated-spectra.gif", fps=None, dpi=None, **kwargs
):
    """
    Create an animation to show how the spectrum changes
    as we flip through every timepoint.

    Parameters
    ----------
    filename : str
        Name of file you'd like to save results in.
        Currently supports only .gif files.
    ax : matplotlib.axes.Axes
        The axes into which this animated plot should go.
    fps : float
        frames/second of animation
    xlim : tuple
        Custom xlimits for the plot
    ylim : tuple
        Custom ylimits for the plot
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap
    scatterkw : dict
        A dictionary of keywords to be passed to `plt.scatter`
    textkw : dict
        A dictionary of keywords to be passed to `plt.text`
    fps : int
        Frames per second for the animation.
    dpi : float, default: :rc:`savefig.dpi`
        Dots per inch for the movie.
    """
    self._setup_animate_spectra(**kwargs)

    # make and save the animation
    animator = ani.FuncAnimation(
        self._animate_spectra_components["fig"],
        self._animate_spectra_components["update"],
        frames=np.arange(0, self.ntime),
        blit=True,
    )
    animator.save(filename, fps=fps, dpi=dpi, savefig_kwargs=dict(facecolor="white"))
    plt.close()


def plot_spectral_resolution(
    self,
    pixels_per_resolution_element=1,
    ax=None,
    w_unit="micron",
    cmap=None,
    vmin=None,
    vmax=None,
    plotkw={},
    **kw,
):
    """
    Plot the spectral resolution as a function of wavelength.

    Parameters
    ----------
    pixels_per_resolution_element : float
        How many pixels do we consider as a resolution element?
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    w_unit : str, astropy.unit.Unit
        The unit for plotting wavelengths.
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength.
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap.
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap.
    plotkw : dict
        A dictionary of keywords passed to `plt.plot`
    """
    R = self.get_spectral_resolution(pixels_per_resolution_element)
    w_unit = u.Unit(w_unit)

    with quantity_support():

        # make sure that the wavelength-based colormap is defined
        self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

        # make sure ax is set up
        if ax is None:
            ax = plt.subplot()
        plt.sca(ax)

        # define default visuals, but let user override with scatterkw
        # this_scatterkw = dict(c=self.get_wavelength_color(self.wavelength))
        # this_scatterkw.update(**scatterkw)

        plt.plot(self.wavelength.to(w_unit), R, **plotkw)
        plt.xlabel(f'Wavelength ({w_unit.to_string("latex_inline")})')
        plt.ylabel(f"$R=\lambda/d\lambda$ ({pixels_per_resolution_element} pixel)")
        return ax

        ax.set_xlabel(f"Time ({self.time.unit.to_string('latex_inline')})")
