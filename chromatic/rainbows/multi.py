from ..imports import *

__all__ = ["MultiRainbow", "compare_rainbows"]


class MultiRainbow:
    """
    MultiRainbow objects collect multiple Rainbow objects together,
    providing a quick interface to apply the same action or visualization
    to all of them at once. It's meant to be a tool to facilitate
    quick comparisons between different pipeline analyses of the
    same dataset.
    """

    def __init__(self, rainbows, names=None):
        """
        Initialize from a list of Rainbows.

        Parameters
        ----------
        rainbows : list
            A list containing two or more Rainbow objects.
        names : list
            A list of names for the Rainbows. These will overwrite
            the names provided in the `.metadata['name']` entry
            of each input Rainbow.
        """

        if not isinstance(rainbows, list):
            raise RuntimeError("Please provide a list of 🌈s.")
        if len(rainbows) < 2:
            raise RuntimeError("Please provide more than one 🌈.")

        # decouple these rainbow from their inputs
        self.rainbows = [x._create_copy() for x in rainbows]

        # keep track of the names
        self.names = names or [(x.name or i) for i, x in enumerate(rainbows)]
        if len(np.unique(self.names)) < len(self.names):
            message = """
            The input names {self.names}
            are not unique. Proceed with caution, or provide
            each 🌈 with a unique name. You can do so either by
            giving each individual 🌈 a name with `rainbow.name='?!?!?!'`
            or by supplying a list of unique names to the `names=`
            keyword argument when starting your comparison.
            """
            cheerfully_suggest(message)

        # make sure the names and rainbows match up
        assert len(self.names) == len(self.rainbows)

        # pull units from first rainbow
        self.w_unit = u.Unit(self.list_of_rainbows[0].wavelength.unit)
        self.t_unit = u.Unit(self.list_of_rainbows[0].time.unit)

    @property
    def list_of_rainbows(self):
        return self.rainbows

    @property
    def dict_of_rainbows(self):
        return {k: v for k, v in zip(self.names, self.rainbows)}

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

    def _setup_panels(self, rows=1, figsize=None, **kw):
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
            constrained_layout=True,
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

        # create a list of interpolating functions for the datasets' dw
        dw_interpolators = [
            interp1d(
                r.wavelength.to_value(self.w_unit),
                np.gradient(r.wavelength.to_value(self.w_unit)),
                fill_value=np.inf,
                bounds_error=False,
            )
            for r in self.list_of_rainbows
        ]

        def smallest_dw(w):
            """
            Helper function to determine the smallest dw
            in the datasets at any particular wavelength.

            Parameters
            ----------
            w : array
                Wavelength (with no units).
            """

            dw = np.inf * np.ones_like(w)
            for interpolator in dw_interpolators:
                dw = np.minimum(dw, interpolator(w))
            return dw

        # compile w and dw, weighted by how many appear in each dataset
        w_represented = np.sort(
            np.hstack(
                [r.wavelength.to_value(self.w_unit) for r in self.list_of_rainbows]
            )
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
            for r in self.list_of_rainbows:
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

        return new_wavelengths * self.w_unit

    def _check_if_wavelengths_are_aligned(self):
        """
        Check whether the wavelength grid is already aligned.

        Returns
        -------
        aligned : bool
            Are the wavelengths the same across all?
        """

        first_rainbow = self.list_of_rainbows[0]

        # check if they're already aligned
        wavelengths_are_same_size = np.all(
            [np.all(r.nwave == first_rainbow.nwave) for r in self.list_of_rainbows]
        )

        if wavelengths_are_same_size:
            wavelengths_are_already_aligned = np.all(
                [
                    np.all(r.wavelength == first_rainbow.wavelength)
                    for r in self.list_of_rainbows
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
            Are the times the same across all Rainbows?
        """
        # check if they're already aligned
        times_are_same_size = np.all(
            [
                np.all(r.ntime == self.list_of_rainbows[0].ntime)
                for r in self.list_of_rainbows
            ]
        )

        if times_are_same_size:
            times_are_already_aligned = np.all(
                [
                    np.all(r.time == self.list_of_rainbows[0].wavelength)
                    for r in self.list_of_rainbows
                ]
            )
        else:
            times_are_already_aligned = False

        return times_are_already_aligned

    @property
    def wavelength(self):
        """
        The shared wavelength axis.
        """
        if self._check_if_wavelengths_are_aligned():
            return self.list_of_rainbows[0].wavelength
        else:
            raise RuntimeError(f"🌈 {self} has more than one wavelength array.")

    @property
    def time(self):
        """
        The shared time axis.
        """
        if self._check_if_times_are_aligned():
            return self.list_of_rainbows[0].time
        else:
            raise RuntimeError(f"🌈 {self} has more than one time array.")

    def align_wavelengths(self, **kw):
        """
        Make a new MultiRainbow where wavelengths have been aligned
        across all the Rainbows, making a guess for a reasonable
        new wavelength grid.

        Parameters
        ----------
        **kw : dict, optional
            Extra keyword arguments will be passed to
            `._guess_good_uniform_wavelength_grid`
        """

        # if the wavelengths are already aligned, don't do anything!
        if self._check_if_wavelengths_are_aligned():
            return self
        else:
            # come up with a decent guess for a shared wavelength grid
            w = self._guess_good_uniform_wavelength_grid(**kw)

            # bin to that shared grid, and don't trim nans at the edges
            return self.bin(wavelength=w, trim=False)

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
        new_rainbows = [r.__getitem__(key) for r in self.list_of_rainbows]
        return MultiRainbow(new_rainbows, names=self.names)

    def normalize(self, **kwargs):
        """
        Normalize by dividing through by the median spectrum and/or lightcurve.

        (see `Rainbow.normalize` for input options)

        Returns
        -------
        normalized : MultiRainbow
            The normalized MultiRainbow.
        """
        new_rainbows = [r.normalize(**kwargs) for r in self.list_of_rainbows]
        return MultiRainbow(new_rainbows, names=self.names)

    def bin(self, **kwargs):
        """
        Bin the rainbow in wavelength and/or time.

        (see `Rainbow.bin` for input options)

        Returns
        -------
        binned : MultiRainbow
            The binned MultiRainbow.
        """
        new_rainbows = [r.bin(**kwargs) for r in self.list_of_rainbows]
        return MultiRainbow(new_rainbows, names=self.names)

    def plot(self, **kwargs):
        """
        Plot flux as sequence of offset light curves.

        (see `Rainbow.plot` for input options)
        """

        # set up a grid of panels
        self._setup_panels()

        # make all the individual plots
        for r, a in zip(self.list_of_rainbows, self.axes):
            r.plot(ax=a, **kwargs)

    def imshow(
        self, vmin=None, vmax=None, rows=1, figsize=None, colorbar=True, **kwargs
    ):
        """
        imshow flux as a function of time (x = time, y = wavelength, color = flux).

        Parameters
        ----------
        vmin : float, optional
            The bottom of the color scale.
        vmax : float, optional
            The top of the color scale.
        rows : int, optional
            Over how many rows should the panels be distributed?
        figsize : tuple, optional
            The (width, height) of the total figure.
        colorbar : bool, optional
            Should we include a colorbar?
        quantity : str, optional
            The fluxlike quantity to imshow.
            (Must be a key of `rainbow.fluxlike`).
        w_unit : str, Unit, optional
            The unit for plotting wavelengths.
        t_unit : str, Unit, optional
            The unit for plotting times.
        aspect : str, optional
            What aspect ratio should be used for the imshow?
        **kwargs : dict, optional
            Most other keyword arguments will be passed on into
            `plt.imshow`. If there's an `plt.imshow` keyword you
            want to set, try it!
        """

        # set up a grid of panels
        self._setup_panels(rows=rows, figsize=figsize)

        # figure out a good shared color limits (unless already supplied)
        vmin = vmin or np.nanmin(
            [np.nanmin(u.Quantity(r.flux).value) for r in self.list_of_rainbows]
        )
        vmax = vmax or np.nanmax(
            [np.nanmax(u.Quantity(r.flux).value) for r in self.list_of_rainbows]
        )

        # make all the individual imshows
        for r, a, i in zip(self.list_of_rainbows, self.axes, range(self.nrainbows)):
            r.imshow(ax=a, vmin=vmin, vmax=vmax, colorbar=False, **kwargs)
            if i == 0:
                xlabel, ylabel = a.get_xlabel(), a.get_ylabel()
            a.set_xlabel("")
            a.set_ylabel("")

        self.figure.supxlabel(xlabel)
        self.figure.supylabel(ylabel)

        # add a colorbar, but steal from all the axes
        if colorbar:
            plt.colorbar(ax=self.axes)

        # set the ylimits to span everything
        plt.ylim(
            np.nanmax([x._imshow_extent[2] for x in self.list_of_rainbows]),
            np.nanmin([x._imshow_extent[3] for x in self.list_of_rainbows]),
        )

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
        filename : str, optional
            Name of file you'd like to save results in.
            Currently supports only .gif files.
        fps : float, optional
            frames/second of animation
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
        textkw : dict, optional
            A dictionary of keywords to be passed to `plt.text`
        """

        # make sure there's the same number of wavelengths for all rainbows
        assert np.all(
            [r.nwave == self.list_of_rainbows[0].nwave for r in self.list_of_rainbows]
        )

        # set up a grid of panels
        self._setup_panels()

        # setup all the individual plots
        for r, a in zip(self.list_of_rainbows, self.axes):
            r._setup_animate_lightcurves(ax=a, **kwargs)

        # the figure should be the same for all panels
        figure = r._animate_lightcurves_components["fi"]

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
            for r in self.list_of_rainbows:
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
        filename : str, optional
            Name of file you'd like to save results in.
            Currently supports only .gif files.
        fps : float, optional
            frames/second of animation
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
        textkw : dict, optional
            A dictionary of keywords to be passed to `plt.text`
        """

        assert np.all(
            [r.ntime == self.list_of_rainbows[0].ntime for r in self.list_of_rainbows]
        )

        # set up a grid of panels
        self._setup_panels()

        # setup all the individual plots
        for r, a in zip(self.list_of_rainbows, self.axes):
            r._setup_animate_spectra(ax=a, **kwargs)

        # the figure should be the same for all panels
        figure = r._animate_spectra_components["fi"]

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
            for r in self.list_of_rainbows:
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


def compare_rainbows(*args, **kwargs):
    """
    Compare a group of Rainbows to each other.

    Parameters
    ----------
    rainbows : list
        A list containing two or more Rainbow objects.

    Returns
    -------
    multirainbow : MultiRainbow
        A MultiRainbow object that contains all of the rainbows to
        be compared. This object can do many of the same actions or
        visualizations as a normal Rainbow, but it will (attempt to)
        automatically apply them to every Rainbow.
    """
    return MultiRainbow(*args, **kwargs)
