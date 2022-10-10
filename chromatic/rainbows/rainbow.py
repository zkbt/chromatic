from ..imports import *
from .readers import *
from ..resampling import *


class Rainbow:
    """
    `Rainbow` (🌈) objects represent brightness as a function
    of both wavelength and time.

    These objects are useful for reading or writing multiwavelength
    time-series datasets in a variety of formats, visualizing these
    data with simple commands, and performing basic calculations.
    `RainbowWithModel` and `SimulatedRainbow` objects inherit from
    `Rainbow`, so all basically all methods and attributes described
    below are available for them too.

    Attributes
    ----------
    wavelike : dict
        A dictionary for quantities with shape `(nwave,),
        for which there's one value for each wavelength.
    timelike : dict
        A dictionary for quantities with shape `(ntime,),
        for which there's one value for each time.
    fluxlike : dict
        A dictionary for quantities with shape `(nwave,ntime),
        for which there's one value for each wavelength and time.
    metadata : dict
        A dictionary containing all other useful information
        that should stay connected to the `Rainbow`, in any format.
    wavelength : Quantity
        The 1D array of wavelengths for this `Rainbow`.
        (This is a property, not an actual attribute.)
    time : Quantity
        The 1D array of times for this `Rainbow`.
        (This is a property, not an actual attribute.)
    flux : array, Quantity
        The 2D array of fluxes for this `Rainbow`.
        (This is a property, not an actual attribute.)
    uncertainty : array, Quantity
        The 2D array of flux uncertainties for this `Rainbow`.
        (This is a property, not an actual attribute.)
    ok : array
        The 2D array of "ok-ness" for this `Rainbow`.
        (This is a property, not an actual attribute.)
    shape : tuple
        The shape of this `Rainbow`'s flux array.
        (This is a property, not an actual attribute.)
    nwave : int
        The number of wavelengths in this `Rainbow`'.
        (This is a property, not an actual attribute.)
    ntime : int
        The number of times in this `Rainbow`'.
        (This is a property, not an actual attribute.)
    nflux : int
        The total number of fluxes in this `Rainbow`' (= `nwave*ntime`).
        (This is a property, not an actual attribute.)
    dt : Quantity
        The typical time offset between adjacent times in this `Rainbow`.
        (This is a property, not an actual attribute.)
    name : array, Quantity
        The name of this `Rainbow`, if one has been set.
        (This is a property, not an actual attribute.)
    """

    # all Rainbows must contain these core dictionaries
    _core_dictionaries = ["fluxlike", "timelike", "wavelike", "metadata"]

    # define which axis is which
    waveaxis = 0
    timeaxis = 1

    # which fluxlike keys will respond to math between objects
    _keys_that_respond_to_math = ["flux"]

    # which keys get uncertainty weighting during binning
    _keys_that_get_uncertainty_weighting = ["flux", "uncertainty"]

    def __init__(
        self,
        filepath=None,
        format=None,
        wavelength=None,
        time=None,
        flux=None,
        uncertainty=None,
        wavelike=None,
        timelike=None,
        fluxlike=None,
        metadata=None,
        name=None,
        **kw,
    ):
        """
        Initialize a `Rainbow` object.

        The `__init__` function is called when a new `Rainbow` is
        instantiated as `r = Rainbow(some, kinds, of=inputs)`.

        The options for inputs are flexible, including the possibility
        to initialize from a file, from arrays with appropriate units,
        from dictionaries with appropriate ingredients, or simply as
        an empty object if no arguments are given.

        Parameters
        ----------
        filepath : str, optional
            The filepath pointing to the file or group of files
            that should be read.
        format : str, optional
            The file format of the file to be read. If None,
            the format will be guessed automatically from the
            filepath.
        wavelength : Quantity, optional
            A 1D array of wavelengths, in any unit.
        time : Quantity, Time, optional
            A 1D array of times, in any unit.
        flux : array, optional
            A 2D array of flux values.
        uncertainty : array, optional
            A 2D array of uncertainties, associated with the flux.
        wavelike : dict, optional
            A dictionary containing 1D arrays with the same
            shape as the wavelength axis. It must at least
            contain the key 'wavelength', which should have
            astropy units of wavelength associated with it.
        timelike : dict, optional
            A dictionary containing 1D arrays with the same
            shape as the time axis. It must at least
            contain the key 'time', which should have
            astropy units of time associated with it.
        fluxlike : dict, optional
            A dictionary containing 2D arrays with the shape
            of (nwave, ntime), like flux. It must at least
            contain the key 'flux'.
        metadata : dict, optional
            A dictionary containing all other metadata
            associated with the dataset, generally lots of
            individual parameters or comments.
        **kw : dict, optional
            Additional keywords will be passed along to
            the function that initializes the rainbow.
            If initializing from arrays (`time=`, `wavelength=`,
            ...), these keywords will be interpreted as
            additional arrays that should be sorted by their
            shape into the appropriate dictionary. If
            initializing from files, the keywords will
            be passed on to the reader.

        Examples
        --------
        Initialize from a file. While this works, a more robust
        solution is probably to use `read_rainbow`, which will
        automatically choose the best of `Rainbow` and `RainbowWithModel`
        ```
        r1 = Rainbow('my-neat-file.abc', format='abcdefgh')
        ```

        Initalize from arrays. The wavelength and time must have
        appropriate units, and the shape of the flux array must
        match the size of the wavelength and time arrays. Other
        arrays that match the shape of any of these quantities
        will be stored in the appropriate location. Other inputs
        not matching any of these will be stored as `metadata.`
        ```
        r2 = Rainbow(
                wavelength=np.linspace(1, 5, 50)*u.micron,
                time=np.linspace(-0.5, 0.5, 100)*u.day,
                flux=np.random.normal(0, 1, (50, 100)),
                some_other_array=np.ones((50,100)),
                some_metadata='wow!'
        )
        ```
        Initialize from dictionaries. The dictionaries must contain
        at least `wavelike['wavelength']`, `timelike['time']`, and
        `fluxlike['flux']`, but any other additional inputs can be
        provided.
        ```
        r3 = Rainbow(
                wavelike=dict(wavelength=np.linspace(1, 5, 50)*u.micron),
                timelike=dict(time=np.linspace(-0.5, 0.5, 100)*u.day),
                fluxlike=dict(flux=np.random.normal(0, 1, (50, 100)))
        )
        ```
        """
        # create a history entry for this action (before other variables are defined)
        h = self._create_history_entry("Rainbow", locals())

        # metadata are arbitrary types of information we need
        self.metadata = {"name": name}

        # wavelike quanities are 1D arrays with nwave elements
        self.wavelike = {}

        # timelike quantities are 1D arrays with ntime elements
        self.timelike = {}

        # fluxlike quantities are 2D arrays with nwave x time elements
        self.fluxlike = {}

        # try to intialize from the exact dictionaries needed
        if (
            (type(wavelike) == dict)
            and (type(timelike) == dict)
            and (type(fluxlike) == dict)
        ):
            self._initialize_from_dictionaries(
                wavelike=wavelike,
                timelike=timelike,
                fluxlike=fluxlike,
                metadata=metadata,
            )
        # then try to initialize from arrays
        elif (wavelength is not None) and (time is not None) and (flux is not None):
            self._initialize_from_arrays(
                wavelength=wavelength,
                time=time,
                flux=flux,
                uncertainty=uncertainty,
                **kw,
            )
            if metadata is not None:
                self.metadata.update(**metadata)
        # then try to initialize from a file
        elif isinstance(filepath, str) or isinstance(filepath, list):
            self._initialize_from_file(filepath=filepath, format=format, **kw)

        # finally, tidy up by guessing the scales
        self._guess_wscale()
        self._guess_tscale()

        # append the history entry to this Rainbow
        self._setup_history()
        self._record_history_entry(h)

    def _sort(self):
        """
        Sort the wavelengths and times, from lowest to highest.
        Attach the unsorted indices to be able to work backwards.
        This sorts the object in-place (not returning a new Rainbow.)

        Returns
        -------
        sorted : Rainbow
            The sorted Rainbow.
        """

        # figure out the indices to sort from low to high
        i_wavelength = np.argsort(self.wavelength)
        i_time = np.argsort(self.time)

        if np.shape(self.flux) != (len(i_wavelength), len(i_time)):
            message = """
            Wavelength, time, and flux arrays don't match;
            the `._sort()` step is being skipped.
            """
            cheerfully_suggest(message)
            return

        if np.any(np.diff(i_wavelength) < 0):
            message = f"""
            The {self.nwave} input wavelengths were not monotonically increasing.
            {self} has been sorted from lowest to highest wavelength.
            If you want to recover the original wavelength order, the original
            wavelength indices are available in `rainbow.original_wave_index`.
            """
            cheerfully_suggest(message)

        if np.any(np.diff(i_time) < 0):
            message = f"""
            The {self.ntime} input times were not monotonically increasing.
            {self} has been sorted from lowest to highest time.
            If you want to recover the original time order, the original
            time indices are available in `rainbow.original_time_index`.
            """
            cheerfully_suggest(message)

        # attach unsorted indices to this array, if the don't exist
        if "original_wave_index" not in self.wavelike:
            self.wavelike["original_wave_index"] = np.arange(self.nwave)
        if "original_time_index" not in self.timelike:
            self.timelike["original_time_index"] = np.arange(self.ntime)

        # sort that copy by wavelength and time
        for k in self.wavelike:
            if self.wavelike[k] is not None:
                self.wavelike[k] = self.wavelike[k][i_wavelength]
        for k in self.timelike:
            if self.timelike[k] is not None:
                self.timelike[k] = self.timelike[k][i_time]
        for k in self.fluxlike:
            if self.fluxlike[k] is not None:
                wave_sorted = self.fluxlike[k][i_wavelength, :]
                self.fluxlike[k][:, :] = wave_sorted[:, i_time]

    def _validate_uncertainties(self):
        """
        Do some checks on the uncertainty values.
        """
        if self.uncertainty is None and len(self.fluxlike) > 0:
            message = f"""
            Hmmm...it's not clear which column corresponds to the
            flux uncertainties for this Rainbow object. The
            available `fluxlike` columns are:
                {self.fluxlike.keys()}
            A long-term solution might be to fix the `from_?!?!?!?`
            reader, but a short-term solution would be to pick one
            of the columns listed above and say something like

            x.fluxlike['uncertainty'] = x.fluxlike['some-other-relevant-error-column']

            where `x` is the Rainbow you just created.
            """
            cheerfully_suggest(message)
            return

        # kludge to replace zero uncertainties
        # if np.all(self.uncertainty == 0):
        #    cheerfully_suggest("\nUncertainties were all 0, replacing them with 1!")
        #        self.fluxlike["uncertainty"] = np.ones_like(self.flux)

    def _initialize_from_dictionaries(
        self, wavelike={}, timelike={}, fluxlike={}, metadata={}
    ):
        """
        Populate from dictionaries in the correct format.

        Parameters
        ----------
        wavelike : dict
            A dictionary containing 1D arrays with the same
            shape as the wavelength axis. It must at least
            contain the key 'wavelength', which should have
            astropy units of wavelength associated with it.
        timelike : dict
            A dictionary containing 1D arrays with the same
            shape as the time axis. It must at least
            contain the key 'time', which should have
            astropy units of time associated with it.
        fluxlike : dict
            A dictionary containing 2D arrays with the shape
            of (nwave, ntime), like flux. It must at least
            contain the key 'flux'.
        metadata : dict
            A dictionary containing all other metadata
            associated with the dataset, generally lots of
            individual parameters or comments.
        """

        # update the three core dictionaries of arrays
        for k in wavelike:
            self.wavelike[k] = wavelike[k] * 1
        for k in timelike:
            self.timelike[k] = timelike[k] * 1
        for k in fluxlike:
            self.fluxlike[k] = fluxlike[k] * 1
        # multiplying by 1 is a kludge to prevent accidental links

        # update the metadata
        self.metadata.update(**metadata)

        # validate that something reasonable got populated
        self._validate_core_dictionaries()

    def _get_core_dictionaries(self):
        """
        Get the core dictionaries of this Rainbow.

        Returns
        -------
        core : dict
            Dictionary containing the keys
            ['wavelike', 'timelike', 'fluxlike', 'metadata']
        """
        return {k: vars(self)[k] for k in self._core_dictionaries}

    def _initialize_from_arrays(
        self, wavelength=None, time=None, flux=None, uncertainty=None, **kw
    ):
        """
        Populate from arrays.

        Parameters
        ----------
        wavelength : Quantity, optional
            A 1D array of wavelengths, in any unit.
        time : Quantity, Time, optional
            A 1D array of times, in any unit.
        flux : array, optional
            A 2D array of flux values.
        uncertainty : array, optional
            A 2D array of uncertainties, associated with the flux.
        **kw : dict, optional
            Additional keywords will be interpreted as arrays
            that should be sorted into the appropriate location
            based on their size.
        """

        # store the wavelength
        self.wavelike["wavelength"] = wavelength * 1

        # store the time
        self.timelike["time"] = time * 1

        # store the flux and uncertainty
        self.fluxlike["flux"] = flux * 1
        if uncertainty is None:
            self.fluxlike["uncertainty"] = np.ones_like(flux) * np.nan
        else:
            self.fluxlike["uncertainty"] = uncertainty * 1

        # sort other arrays by shape
        for k, v in kw.items():
            self._put_array_in_right_dictionary(k, v)

        # validate that something reasonable got populated
        self._validate_core_dictionaries()

    def _put_array_in_right_dictionary(self, k, v):
        """
        Sort an input into the right core dictionary
        (timelike, wavelike, fluxlike) based on its shape.

        Parameters
        ----------
        k : str
            The key for the (appropriate) dictionary.
        v : array
            The quantity to sort.
        """
        if np.shape(v) == self.shape:
            self.fluxlike[k] = v * 1
        elif np.shape(v) == (self.nwave,):
            self.wavelike[k] = v * 1
        elif np.shape(v) == (self.ntime,):
            self.timelike[k] = v * 1
        else:
            raise ValueError(f"'{k}' doesn't fit anywhere!")

    def _initialize_from_file(self, filepath=None, format=None, **kw):
        """
        Populate from a filename or group of files.

        Parameters
        ----------
        filepath : str, optional
            The filepath pointing to the file or group of files
            that should be read.
        format : str, optional
            The file format of the file to be read. If None,
            the format will be guessed automatically from the
            filepath.
        **kw : dict, optional
            Additional keywords will be passed on to the reader.
        """

        # make sure we're dealing with a real filename
        assert filepath is not None

        # pick the appropriate reader
        reader = guess_reader(filepath=filepath, format=format)
        reader(self, filepath, **kw)

        # validate that something reasonable got populated
        self._validate_core_dictionaries()
        self._validate_uncertainties()
        self._guess_wscale()
        self._guess_tscale()

    def _create_copy(self):
        """
        Create a copy of self, with the core dictionaries copied.
        """
        new = type(self)()
        new._initialize_from_dictionaries(
            **copy.deepcopy(self._get_core_dictionaries())
        )
        return new

    def _guess_wscale(self, relative_tolerance=0.05):
        """
        Try to guess the wscale from the wavelengths.

        Parameters
        ----------

        relative_tolerance : float
            The fractional difference to which the differences
            between wavelengths should match in order for a
            linear or logarithmic wavelength scale to be
            assigned. For example, the default value of 0.01
            means that the differences between all wavelength
            must be within 1% of each other for the wavelength
            scale to be called linear.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # give up if there's no wavelength array
            if self.wavelength is None:
                return "?"

            # calculate difference arrays
            w = self.wavelength.value
            dw = np.diff(w)
            dlogw = np.diff(np.log(w))

            # test the three options
            if np.allclose(dw, np.median(dw), rtol=relative_tolerance):
                self.metadata["wscale"] = "linear"
            elif np.allclose(dlogw, np.median(dlogw), rtol=relative_tolerance):
                self.metadata["wscale"] = "log"
            else:
                self.metadata["wscale"] = "?"

    def _guess_tscale(self, relative_tolerance=0.05):
        """
        Try to guess the tscale from the times.

        Parameters
        ----------

        relative_tolerance : float
            The fractional difference to which the differences
            between times should match in order for us to call
            the times effectively uniform, or for us to treat
            them more carefully as an irregular or gappy grid.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # give up if there's no time array
            if self.time is None:
                return "?"

            # calculate difference arrays
            t = self.time.value
            dt = np.diff(t)
            with warnings.catch_warnings():
                # (don't complain about negative time)
                warnings.simplefilter("ignore")
                dlogt = np.diff(np.log(t))

            # test the three options
            if np.allclose(dt, np.median(dt), rtol=relative_tolerance):
                self.metadata["tscale"] = "linear"
            # elif np.allclose(dlogt, np.median(dlogt), rtol=relative_tolerance):
            #    self.metadata["tscale"] = "log"
            else:
                self.metadata["tscale"] = "?"

    @property
    def name(self):
        """
        The name of this `Rainbow` object.
        """
        return self.metadata.get("name", None)

    @property
    def wavelength(self):
        """
        The 1D array of wavelengths (with astropy units of length).
        """
        return self.wavelike.get("wavelength", None)

    @property
    def time(self):
        """
        The 1D array of time (with astropy units of time).
        """
        return self.timelike.get("time", None)

    @property
    def flux(self):
        """
        The 2D array of fluxes (row = wavelength, col = time).
        """
        return self.fluxlike.get("flux", None)

    @property
    def uncertainty(self):
        """
        The 2D array of uncertainties on the fluxes.
        """
        return self.fluxlike.get("uncertainty", None)

    @property
    def ok(self):
        """
        The 2D array of whether data is OK (row = wavelength, col = time).
        """

        # assemble from three possible arrays
        ok = self.fluxlike.get("ok", np.ones(self.shape).astype(bool))
        ok = (
            ok
            * self.wavelike.get("ok", np.ones(self.nwave).astype(bool))[:, np.newaxis]
        )
        ok = (
            ok
            * self.timelike.get("ok", np.ones(self.ntime).astype(bool))[np.newaxis, :]
        )

        # make sure flux is finite
        if self.flux is not None:
            ok = ok * np.isfinite(self.flux)

        # weird kludge to deal with rounding errors (particularly in two-step .bin)
        if ok.dtype == bool:
            return ok
        elif np.all((ok == 1) | (ok == 0)):
            return ok.astype(bool)
        else:
            return np.round(ok, decimals=12)

    @property
    def _time_label(self):
        return self.metadata.get("time_label", "Time")

    @property
    def _wave_label(self):
        return self.metadata.get("wave_label", "Wavelength")

    def __getattr__(self, key):
        """
        If an attribute/method isn't explicitly defined,
        try to pull it from one of the core dictionaries.

        Let's say you want to get the 2D uncertainty array
        but don't want to type `self.fluxlike['uncertainty']`.
        You could instead type `self.uncertainty`, and this
        would try to search through the four standard
        dictionaries to pull out the first `uncertainty`
        it finds.

        Parameters
        ----------
        key : str
            The attribute we're trying to get.
        """
        if key not in self._core_dictionaries:
            for dictionary_name in self._core_dictionaries:
                try:
                    return self.__dict__[dictionary_name][key]
                except KeyError:
                    pass
        message = f"🌈.{key} does not exist for this Rainbow"
        raise AttributeError(message)

    def __setattr__(self, key, value):
        """
        When setting a new attribute, try to sort it into the
        appropriate core directory based on its size.

        Let's say you have some quantity that has the same
        shape as the wavelength array and you'd like to attach
        it to this Rainbow object. This will try to save it
        in the most relevant core dictionary (of the choices
        timelike, wavelike, fluxlike).

        Parameters
        ----------
        key : str
            The attribute we're trying to get.
        value : array
            The quantity we're trying to attach to that name.
        """
        try:
            if key in self._core_dictionaries:
                raise ValueError("Trying to set a core dictionary.")
            elif key == "wavelength":
                self.wavelike["wavelength"] = value * 1
                self._validate_core_dictionaries()
            elif key == "time":
                self.timelike["time"] = value * 1
                self._validate_core_dictionaries()
            elif key in ["flux", "uncertainty", "ok"]:
                self.fluxlike[key] = value * 1
                self._validate_core_dictionaries()
            elif isinstance(value, str):
                self.metadata[key] = value
            else:
                self._put_array_in_right_dictionary(key, value)
        except (AttributeError, ValueError):
            self.__dict__[key] = value

    @property
    def _nametag(self):
        """
        This short phrase will preface everything
        said with `self.speak()`.
        """
        return f"🌈({self.nwave}w, {self.ntime}t)"

    @property
    def shape(self):
        """
        The shape of the flux array (nwave, ntime).
        """
        return (self.nwave, self.ntime)

    @property
    def nwave(self):
        """
        The number of wavelengths.
        """
        if self.wavelength is None:
            return 0
        else:
            return len(self.wavelength)

    @property
    def ntime(self):
        """
        The number of times.
        """
        if self.time is None:
            return 0
        else:
            return len(self.time)

    @property
    def dt(self):
        """
        The typical timestep.
        """
        if self.time is None:
            return None
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return np.nanmedian(np.diff(self.time)).to(u.minute)

    @property
    def nflux(self):
        """
        The total number of fluxes.
        """
        return np.prod(self.shape)

    def _validate_core_dictionaries(self):
        """
        Do some simple checks to make sure this Rainbow
        is populated with the minimal data needed to do anything.
        It shouldn't be run before the Rainbow is fully
        initialized; otherwise, it might complain about
        a half-populated object.
        """

        # make sure there are some times + wavelengths defined
        if self.ntime is None:
            cheerfully_suggest(
                f"""
            No times are defined for this Rainbow.
            """
            )
        if self.nwave is None:
            cheerfully_suggest(
                f"""
            No wavelengths are defined for this Rainbow.
            """
            )

        # warn if the times and wavelengths are the same size
        if (self.nwave == self.ntime) and (self.ntime is not None) and (self.ntime > 1):
            cheerfully_suggest(
                f"""
            The number of wavelengths ({self.nwave}) is the same as the
            number of times ({self.ntime}). This is fine, we suppose
            (<mock exasperated sigh>), but here are few reasons you might
            want to reconsider letting them have the same size:
                (1) Mathemetical operations and variabile assignment
                    inside this Rainbow make guesses about whether a quantity
                    is wavelike or timelike based on its shape; these features
                    will fail (or even worse do something mysterious) if
                    there are the same numbers of wavelengths and times.
                (2) For your own darn sake, if your fluxlike arrays are
                    all square, it's going to be very easy for you to accidentally
                    transpose them and not realize it.
                (3) It's very unlikely that your real data had exactly the same
                    number of times and wavelengths, so we're guessing that you
                    probably just created these arrays from scratch, which
                    hopefully means it's not too annoying to just make them
                    have different numbers of wavelengths and times.
            Thanks!
            """
            )

        # does the flux have the right shape?
        if self.shape != np.shape(self.flux):
            message = f"""
            Something doesn't line up!
            The flux array has a shape of {np.shape(self.flux)}.
            The wavelength array has {self.nwave} wavelengths.
            The time array has {self.ntime} times.
            """
            if self.shape == np.shape(self.flux)[::-1]:
                cheerfully_suggest(
                    f"""{message}
                    Any chance your flux array is transposed?
                    """
                )
            else:
                cheerfully_suggest(message)

        for n in ["uncertainty", "ok"]:
            x = getattr(self, n)
            if x is not None:
                if x.shape != np.shape(self.flux):
                    message = f"""
                    Watch out! The '{n}' array has
                    a shape of {x.shape}, which doesn't match the
                    flux array's shape of {np.shape(self.flux)}.
                    """
                    cheerfully_suggest(message)

        # make sure 2D arrays are uniquely named from 1D
        for k in tuple(self.fluxlike.keys()):
            if (k in self.wavelike) or (k in self.timelike):
                self.fluxlike[f"{k}_2d"] = self.fluxlike.pop(k)

        if "ok" in self.fluxlike:
            is_nan = np.isnan(self.fluxlike["flux"])
            self.fluxlike["ok"][is_nan] = 0

        # make sure no arrays are accidentally pointed to each other
        # (if they are, sorting will get really messed up!)
        for d in ["fluxlike", "wavelike", "timelike"]:
            core_dictionary = self.get(d)
            for k1, v1 in core_dictionary.items():
                for k2, v2 in core_dictionary.items():
                    if k1 != k2:
                        assert v1 is not v2

        self._sort()

    def _make_sure_wavelength_edges_are_defined(self):
        """
        Make sure there are some wavelength edges defined.
        """
        if self.nwave <= 1:
            return
        if ("wavelength_lower" not in self.wavelike) or (
            "wavelength_upper" not in self.wavelike
        ):
            if self.metadata.get("wscale", None) == "log":
                l, u = calculate_bin_leftright(np.log(self.wavelength.value))
                self.wavelike["wavelength_lower"] = np.exp(l) * self.wavelength.unit
                self.wavelike["wavelength_upper"] = np.exp(u) * self.wavelength.unit
            elif self.metadata.get("wscale", None) == "linear":
                l, u = calculate_bin_leftright(self.wavelength)
                self.wavelike["wavelength_lower"] = l
                self.wavelike["wavelength_upper"] = u
            else:
                l, u = calculate_bin_leftright(self.wavelength)
                self.wavelike["wavelength_lower"] = l
                self.wavelike["wavelength_upper"] = u

    def _make_sure_time_edges_are_defined(self, redo=True):
        """
        Make sure there are some time edges defined.
        """
        if self.ntime <= 1:
            return
        if (
            ("time_lower" not in self.timelike)
            or ("time_upper" not in self.timelike)
            or redo
        ):
            if self.metadata.get("tscale", None) == "log":
                lower, upper = calculate_bin_leftright(np.log(self.time.value))
                self.timelike["time_lower"] = np.exp(lower) * self.time.unit
                self.timelike["time_upper"] = np.exp(upper) * self.time.unit
            else:
                lower, upper = calculate_bin_leftright(self.time)
                self.timelike["time_lower"] = lower
                self.timelike["time_upper"] = upper

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

        i_wavelength, i_time = key
        # create a history entry for this action (before other variables are defined)
        h = self._create_history_entry("__getitem__", locals())

        # create a copy
        new = self._create_copy()

        # make sure we don't drop down to 1D arrays
        if isinstance(i_wavelength, int):
            i_wavelength = [i_wavelength]

        if isinstance(i_time, int):
            i_time = [i_time]

        # do indexing of wavelike
        for w in self.wavelike:
            new.wavelike[w] = self.wavelike[w][i_wavelength]

        # do indexing of timelike
        for t in self.timelike:
            new.timelike[t] = self.timelike[t][i_time]

        # do indexing of fluxlike
        for f in self.fluxlike:
            # (indexing step by step seems more stable)
            if self.fluxlike[f] is None:
                continue
            temporary = self.fluxlike[f][i_wavelength, :]
            new.fluxlike[f] = temporary[:, i_time]

        # finalize the new rainbow
        new._validate_core_dictionaries()
        new._guess_wscale()
        new._guess_tscale()

        # append the history entry to the new Rainbow
        new._record_history_entry(h)

        return new

    def __repr__(self):
        """
        How should this object be represented as a string?
        """
        n = self.__class__.__name__.replace("Rainbow", "🌈")
        if self.name is not None:
            n += f"'{self.name}'"
        return f"<{n}({self.nwave}w, {self.ntime}t)>"

    # import the basic operations for Rainbows
    from .actions.operations import (
        _apply_operation,
        _broadcast_to_fluxlike,
        _raise_ambiguous_shape_error,
        __add__,
        __sub__,
        __mul__,
        __truediv__,
        __eq__,
        diff,
    )

    # import other actions that return other Rainbows
    from .actions import (
        normalize,
        _is_probably_normalized,
        bin,
        bin_in_time,
        bin_in_wavelength,
        trim,
        trim_times,
        trim_wavelengths,
        shift,
        _create_shared_wavelength_axis,
        align_wavelengths,
        inject_transit,
        inject_systematics,
        inject_noise,
        inject_spectrum,
        inject_outliers,
        flag_outliers,
        fold,
        mask_transit,
        compare,
        get_average_lightcurve_as_rainbow,
        get_average_spectrum_as_rainbow,
        _create_fake_wavelike_quantity,
        _create_fake_timelike_quantity,
        _create_fake_fluxlike_quantity,
        remove_trends,
        attach_model,
        inflate_uncertainty,
    )

    # import summary statistics for each wavelength
    from .get.wavelike import (
        get_average_spectrum,
        get_median_spectrum,
        get_spectral_resolution,
        get_expected_uncertainty,
        get_measured_scatter,
        get_measured_scatter_in_bins,
        get_for_wavelength,
        get_ok_data_for_wavelength,
    )

    # import summary statistics for each wavelength
    from .get.timelike import (
        get_average_lightcurve,
        get_median_lightcurve,
        get_for_time,
        get_ok_data_for_time,
        get_times_as_astropy,
        set_times_from_astropy,
    )

    # import visualizations that can act on Rainbows
    from .visualizations import (
        imshow,
        pcolormesh,
        plot_lightcurves,
        _setup_animate_lightcurves,
        animate_lightcurves,
        _setup_animate_spectra,
        animate_spectra,
        _setup_animated_scatter,
        setup_wavelength_colors,
        _make_sure_cmap_is_defined,
        get_wavelength_color,
        imshow_quantities,
        plot_quantities,
        imshow_interact,
        plot_spectra,
        plot,
        plot_histogram,
        _scatter_timelike_or_wavelike,
        _get_plot_directory,
        _label_plot_file,
        savefig,
    )

    from .visualizations.wavelike import (
        plot_spectral_resolution,
        plot_noise_comparison,
        plot_noise_comparison_in_bins,
        plot_average_spectrum,
        plot_median_spectrum,
    )

    from .visualizations.timelike import plot_average_lightcurve, plot_median_lightcurve

    from .converters import (
        to_nparray,
        to_df,
    )

    from .helpers import (
        _setup_history,
        _record_history_entry,
        _remove_last_history_entry,
        _create_history_entry,
        history,
        help,
        save,
        get,
    )
