from ..imports import *


class MultiRainbow:
    def __init__(self, list_of_rainbows):
        self.rainbows = list_of_rainbows

    def __repr__(self):
        return f"<{self.rainbows}>"

    @property
    def nrainbows(self):
        return len(self.rainbows)

    def _setup_panels(self, figsize=None):
        if figsize == None:
            default_figsize = plt.matplotlib.rcParams["figure.figsize"]
            figsize = [self.nrainbows * default_figsize[0], default_figsize[1]]
        self.figure, self.axes = plt.subplots(
            1, self.nrainbows, sharex=True, sharey=True, figsize=figsize
        )

    def normalize(self, **kwargs):
        """
        Normalize by dividing through by the median spectrum and/or lightcurve.

        Parameters
        ----------
        wavelength : bool
            Should we divide by the median spectrum?

        time : bool
            Should we divide by the median light curve?

        Returns
        -------
        normalized : MultiRainbow
            The normalized MultiRainbow.
        """
        new_rainbows = [r.normalize(**kwargs) for r in self.rainbows]
        return MultiRainbow(new_rainbows)

    def bin(self, **kwargs):
        """
        Bin the rainbow in wavelength and/or time.

        The time-setting order of precendence is
        [`time`, `dt`], meaning that if `time` is set,
        any values given for `dt` will be ignored.
        The wavelength-setting order of precendence is
        [`wavelength`, `dw`, `R`], meaning that if `wavelength`
        is set any values of `dw` or `R` will be ignored, and
        if `dw` is set any value of `R` will be ignored.

        Parameters
        ----------
        dt : astropy.units.Quantity
            The d(time) bin size for creating a grid
            that is uniform in linear space.
        time : array of astropy.units.Quantity
            An array of times, if you just want to give
            it an entirely custom array.
        R : float
            The spectral resolution for creating a grid
            that is uniform in logarithmic space.
        dw : astropy.units.Quantity
            The d(wavelength) bin size for creating a grid
            that is uniform in linear space.
        wavelength : array of astropy.units.Quantity
            An array of wavelengths, if you just want to give
            it an entirely custom array.

        Returns
        -------
        binned : MultiRainbow
            The binned MultiRainbow.
        """
        new_rainbows = [r.bin(**kwargs) for r in self.rainbows]
        return MultiRainbow(new_rainbows)

    def plot(self, **kwargs):
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

        # set up a grid of panels
        self._setup_panels()

        # make all the individual plots
        for r, a in zip(self.rainbows, self.axes):
            r.plot(ax=a, **kwargs)

    def imshow(self, vmin=None, vmax=None, **kwargs):
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

        # set up a grid of panels
        self._setup_panels()

        # figure out a good shared color limits
        vmin = vmin or np.min([np.nanmin(r.flux) for r in self.rainbows])
        vmax = vmax or np.max([np.nanmax(r.flux) for r in self.rainbows])

        # make all the individual imshows
        for r, a in zip(self.rainbows, self.axes):
            r.imshow(ax=a, vmin=vmin, vmax=vmax, **kwargs)

    def animate_lightcurves(
        self,
        filename="multi-animated-lightcurves.gif",
        fps=None,
        dpi=None,
        **kwargs,
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

        assert np.all([r.nwave == self.rainbows[0].nwave for r in self.rainbows])

        # set up a grid of panels
        self._setup_panels()

        # setup all the individual plots
        for r, a in zip(self.rainbows, self.axes):
            r._setup_animate_lightcurves(ax=a, **kwargs)

        # the figure should be the same for all panels
        figure = r._animate_lightcurves_components["fig"]

        def update(frame):
            """
            This function will be called to update each frame
            of the animation.

            Parameters
            ----------
            frame : int
                An integer that will advance with each frame.
            """
            artists = []
            for r in self.rainbows:
                r._animate_lightcurves_components["update"](frame)
                for k in ["text", "scatter"]:
                    artists.append(r._animate_lightcurves_components[k])
            return artists

        # make and save the animation
        animator = ani.FuncAnimation(
            figure,
            update,
            frames=np.arange(0, r.nwave),
            blit=True,
        )
        animator.save(
            filename, fps=fps, dpi=dpi, savefig_kwargs=dict(facecolor="white")
        )
        plt.close()

    def animate_spectra(
        self, filename="multi-animated-spectra.gif", fps=None, dpi=None, **kwargs
    ):
        """
        Create an animation to show how the spectrum changes
        as we flip through every timepoint.

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

        assert np.all([r.ntime == self.rainbows[0].ntime for r in self.rainbows])

        # set up a grid of panels
        self._setup_panels()

        # setup all the individual plots
        for r, a in zip(self.rainbows, self.axes):
            r._setup_animate_spectra(ax=a, **kwargs)

        # the figure should be the same for all panels
        figure = r._animate_spectra_components["fig"]

        def update(frame):
            """
            This function will be called to update each frame
            of the animation.

            Parameters
            ----------
            frame : int
                An integer that will advance with each frame.
            """
            artists = []
            for r in self.rainbows:
                r._animate_spectra_components["update"](frame)
                for k in ["text", "scatter"]:
                    artists.append(r._animate_spectra_components[k])
            return artists

        # make and save the animation
        animator = ani.FuncAnimation(
            figure,
            update,
            frames=np.arange(0, r.ntime),
            blit=True,
        )
        animator.save(
            filename, fps=fps, dpi=dpi, savefig_kwargs=dict(facecolor="white")
        )
        plt.close()
