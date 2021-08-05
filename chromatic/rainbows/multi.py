from ..imports import *


class MultiRainbow:
    """
    MultiRainbow objects collect multiple Rainbow objects together,
    providing a quick interface to apply the same action or visualization
    to all of them at once. It's meant to be a tool to facilitate
    quick comparisons between different pipeline analyses of the
    same dataset.
    """

    def __init__(self, list_of_rainbows, names=None):
        """
        Initialize from a list of Rainbows.

        Parameters
        ----------
        list_of_rainbows : list
            A list containing 1 or more Rainbow objects.
        names : list
            A list of names with which to label the Rainbows.
        """

        # store the Rainbow objects
        self.rainbows = list_of_rainbows

        # the list of names associated with those objects
        self.names = names

        # make sure the names and rainbows match up
        if self.names is not None:
            assert len(self.names) == len(self.rainbows)

    def __repr__(self):
        """
        How should this object be represented?
        """
        return f"<MultiRainbow({self.rainbows})>"

    @property
    def nrainbows(self):
        """
        What's the number of rainbows in here?
        """
        return len(self.rainbows)

    def _setup_panels(self, rows=1, figsize=None):
        """
        Set up a grid of panels to plot into.

        Parameters
        ----------
        rows : int
            The number of rows into which to arrange the
            plot panels. (This will probably only be needed
            for working with more than about 3 Rainbows.)
        figsize : tuple
            The size (width, height) for the figure to create.
            If left blank, this will be estimated from the
            default matplotlib figsize and the number of
            figures being produced.
        """

        # estimate a figure size from current matplotlib defaults
        if figsize == None:
            default_figsize = plt.matplotlib.rcParams["figure.figsize"]
            figsize = [self.nrainbows * default_figsize[0], default_figsize[1]]

        # create the figure and grid of axes as subplots
        self.figure, self.axes = plt.subplots(
            rows,
            int(np.ceil(self.nrainbows / rows)),
            sharex=True,
            sharey=True,
            figsize=figsize,
        )

        # give titles to the plot panels
        if self.names is not None:
            for name, ax in zip(self.names, self.axes):
                ax.set_title(name)

    def _guess_good_uniform_wavelength_grid(
        self, fraction_to_supersample=0.5, wscale=None, plot=False
    ):
        """
        Guess a good shared wavelength grid that would work
        reasonably well for all the Rainbows in this set.

        Parameters
        ----------
        fraction_to_supersample : float
            When setting the new wavelength bins, what fraction of bins
            in the original Rainbows will be super-sampled? 0.0 will create
            new bins that are just larger than the largest original bin,
            1.0 will create new bins that are just smaller than the smallest
            original bin, and 0.5 is the default.

        plot : bool
            Should we make a plot of the d[wavelength]/d[bin], for diagnostics?
        """

        # pick a unit to make sure everything's consistent
        w_unit = u.micron

        # create a list of interpolating functions for the datasets' dw
        dw_interpolators = [
            interp1d(
                r.wavelength.to(w_unit).value,
                np.gradient(r.wavelength.to(w_unit).value),
                fill_value=np.inf,
                bounds_error=False,
            )
            for r in self.rainbows
        ]

        def smallest_dw(w):
            """
            Helper function to determine the smallest dw
            in the datasets at any particular wavelength.

            Parameters
            ----------
            w : np.array
                Wavelength (with no units).
            """

            dw = np.inf * np.ones_like(w)
            for interpolator in dw_interpolators:
                dw = np.minimum(dw, interpolator(w))
            return dw

        # compile w and dw, weighted by how many appear in each dataset
        w_represented = np.sort(
            np.hstack([r.wavelength.to(w_unit).value for r in self.rainbows])
        )
        dw_represented = smallest_dw(w_represented)

        # fit a line to log(dw) as a function of log(w)
        slope, intercept = np.polyfit(np.log(w_represented), np.log(dw_represented), 1)
        dw_model = np.exp(np.polyval([slope, intercept], np.log(w_represented)))

        # decide if it should be linear or logarithmic
        if wscale is None:
            if slope < 0.5:
                wscale = "linear"
            else:
                wscale = "log"

        if wscale == "linear":
            dw = np.percentile(dw_represented, 100 * (1 - fraction_to_supersample))
            new_wavelengths = np.arange(
                np.nanmin(w_represented), np.nanmax(w_represented), dw
            )
        elif wscale == "log":
            R = np.percentile(
                w_represented / dw_represented, 100 * fraction_to_supersample
            )
            new_wavelengths = np.exp(
                np.arange(
                    np.log(np.nanmin(w_represented)),
                    np.log(np.nanmax(w_represented)),
                    1 / R,
                )
            )
        else:
            raise RuntimeError('🌈 Currently only "log" or "linear" are supported.')

        if plot:

            # plot the individual datasets
            for r in self.rainbows:
                plt.plot(
                    r.wavelength, np.gradient(r.wavelength), marker="o", alpha=0.25
                )

            # plot the minimum dw at each wavelength
            plt.plot(
                w_represented, dw_represented, marker="o", color="black", alpha=0.25
            )

            # convert to log scale, so R=constant is a slope=1 line
            plt.yscale("log")
            plt.xscale("log")
            plt.axis("scaled")

            # plot the line fit result
            plt.plot(
                w_represented,
                dw_model,
                linestyle="--",
                color="black",
                linewidth=2,
                alpha=0.25,
            )

            # plot the new logarithmic or linear wavelength grid
            plt.plot(
                new_wavelengths,
                np.gradient(new_wavelengths),
                color="black",
                linewidth=20,
                alpha=0.25,
            )
            plt.title(f"{slope:.3} -> {wscale}!")
            plt.ylabel("d[wavelength]/d[bin]")

        return new_wavelengths * w_unit

    def _check_if_wavelengths_are_aligned(self):
        """
        Check whether the wavelength grid is already aligned.

        Returns
        -------
        aligned : bool
            Are the wavelengths the same across all?
        """
        # check if they're already aligned
        wavelengths_are_same_size = np.all(
            [np.all(r.nwave == self.rainbows[0].nwave) for r in self.rainbows]
        )

        if wavelengths_are_same_size:
            wavelengths_are_already_aligned = np.all(
                [
                    np.all(r.wavelength == self.rainbows[0].wavelength)
                    for r in self.rainbows
                ]
            )
        else:
            wavelengths_are_already_aligned = False

        return wavelengths_are_already_aligned

    def _check_if_times_are_aligned(self):
        """
        Check whether the time grid is already aligned.

        Returns
        -------
        aligned : bool
            Are the times the same across all?
        """
        # check if they're already aligned
        times_are_same_size = np.all(
            [np.all(r.ntime == self.rainbows[0].ntime) for r in self.rainbows]
        )

        if times_are_same_size:
            times_are_already_aligned = np.all(
                [np.all(r.time == self.rainbows[0].wavelength) for r in self.rainbows]
            )
        else:
            times_are_already_aligned = False

        return times_are_already_aligned

    @property
    def wavelength(self):
        if self._check_if_wavelengths_are_aligned():
            return self.rainbows[0].wavelength
        else:
            raise RuntimeError(f"🌈 {self} has more than one wavelength array.")

    @property
    def time(self):
        if self._check_if_times_are_aligned():
            return self.rainbows[0].time
        else:
            raise RuntimeError(f"🌈 {self} has more than one time array.")

    def align_wavelengths(self, **kw):
        """
        Make a new MultiRainbow where wavelengths have
        been aligned across all the Rainbows, making
        a guesst for whether
        """

        if self._check_if_wavelengths_are_aligned():
            return self
        else:
            w = self._guess_good_uniform_wavelength_grid(**kw)
            return self.bin(wavelength=w)

    def __getitem__(self, key):
        """
        Trim a rainbow by indexing, slicing, or masking.
        Two indices must be provided (`[:,:]`).

        Examples
        --------
        ```
        r[:,:]
        r[10:20, :]
        r[np.arange(10,20), :]
        r[r.wavelength > 1*u.micron, :]
        r[:, np.abs(r.time) < 1*u.hour]
        r[r.wavelength > 1*u.micron, np.abs(r.time) < 1*u.hour]
        ```

        Parameters
        ----------
        key : tuple
            The (wavelength, time) slices, indices, or masks.
        """
        new_rainbows = [r.__getitem__(key) for r in self.rainbows]
        return MultiRainbow(new_rainbows, names=self.names)

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
        return MultiRainbow(new_rainbows, names=self.names)

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
        return MultiRainbow(new_rainbows, names=self.names)

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
        vmin = vmin or np.nanmin([np.nanmin(r.flux) for r in self.rainbows])
        vmax = vmax or np.nanmax([np.nanmax(r.flux) for r in self.rainbows])

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
