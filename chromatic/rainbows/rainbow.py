from ..imports import *
from ..talker import Talker
from ..resampling import *
from .readers import *


class Rainbow(Talker):
    """
    Rainbow objects represent the flux of some astrophysical
    object as a function of both wavelength and time.

    The general stucture of a Rainbow object contains a
    1D array of wavelengths, a 1D array of times, and a
    2D array of flux values (with row as wavelength and
    column as time).
    """

    # all Rainbows must contain these core dictionaries
    _core_dictionaries = ["fluxlike", "timelike", "wavelike", "metadata"]

    # define which axis is which
    waveaxis = 0
    timeaxis = 1

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
        **kw,
    ):
        """
        Initialize a Rainbow object from a file, from arrays
        with appropriate units, or from dictionaries with
        appropriate ingredients.

        Parameters
        ----------

        filepath : str
            The filepath pointing to the file or group of files
            that should be read.
        format : str
            The file format of the file to be read. If None,
            the format will be guessed automatically from the
            filepath.

        wavelength : astropy.unit.Quantity
            A 1D array of wavelengths, in any unit.
        time : astropy.unit.Quantity or astropy.unit.Time
            A 1D array of times, in any unit.
        flux : np.array
            A 2D array of flux values.
        uncertainty : np.array
            A 2D array of uncertainties, associated with the flux.

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

        Examples
        --------
        ```
        # initalize from a file
        r1 = Rainbow('my-neat-file.abc', format='abcdefgh')

        # initalize from arrays
        r2 = Rainbow(wavelength=np.linspace(1, 5, 50)*u.micron,
                     time=np.linspace(-0.5, 0.5, 100)*u.day,
                     flux=np.random.normal(0, 1, (50, 100)))

        # initialize from dictionaries
        f3 = Rainbow(wavelike=dict(wavelength=np.linspace(1, 5, 50)*u.micron),
                     timelike=dict(time=np.linspace(-0.5, 0.5, 100)*u.day),
                     fluxlike=dict(flux=np.random.normal(0, 1, (50, 100))))

        """
        # metadata are arbitrary types of information we need
        self.metadata = {}

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
        # then try to initialize from a file
        elif (type(filepath) == str) or (type(filepath) == list):
            self._initialize_from_file(filepath=filepath, format=format)

        # finally, tidy up by guessing the wavelength scale
        self._guess_wscale()

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

        # update the four core dictionaries
        self.wavelike.update(**wavelike)
        self.timelike.update(**timelike)
        self.fluxlike.update(**fluxlike)
        self.metadata.update(**metadata)

        # validate that something reasonable got populated
        self._validate_core_dictionaries()

    def _initialize_from_arrays(
        self, wavelength=None, time=None, flux=None, uncertainty=None, **kw
    ):
        """
        Populate from arrays.

        Parameters
        ----------
        wavelength : astropy.unit.Quantity
            A 1D array of wavelengths, in any unit.
        time : astropy.unit.Quantity or astropy.unit.Time
            A 1D array of times, in any unit.
        flux : np.array
            A 2D array of flux values.
        uncertainty : np.array
            A 2D array of uncertainties, associated with the flux.
        """

        # store the wavelength
        self.wavelike["wavelength"] = wavelength

        # store the time
        self.timelike["time"] = time

        # store the flux and uncertainty
        self.fluxlike["flux"] = flux
        self.fluxlike["uncertainty"] = uncertainty

        # validate that something reasonable got populated
        self._validate_core_dictionaries()

    def _initialize_from_file(self, filepath=None, format=None):
        """
        Populate from a filename or group of files.

        Parameters
        ----------

        filepath : str
            The filepath pointing to the file or group of files
            that should be read.
        format : str
            The file format of the file to be read. If None,
            the format will be guessed automatically from the
            filepath.
        """

        # make sure we're dealing with a real filename
        assert filepath is not None

        # pick the appropriate reader
        reader = guess_reader(filepath=filepath, format=format)
        reader(self, filepath)

        # validate that something reasonable got populated
        self._validate_core_dictionaries()

    def _guess_wscale(self, relative_tolerance=0.01):
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

    def _guess_tscale(self, relative_tolerance=0.01):
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

        # give up if there's no time array
        if self.time is None:
            return "?"

        # calculate difference arrays
        t = self.time.value
        dt = np.diff(t)

        # test the three options
        if np.allclose(dt, np.median(dt), rtol=relative_tolerance):
            self.metadata["tscale"] = "uniform"
        else:
            self.metadata["wscale"] = "?"

    @property
    def wavelength(self):
        """
        The 1D array of wavelengths (with astropy units of length).
        """
        return self.wavelike.get("wavelength", None)

    @property
    def time(self):
        """
        The 1D array of wavelengths (with astropy units of length).
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
        raise AttributeError()

    # TODO - what should we do with __setattr__?
    #   actually allow to reset things in metadata?
    #   give a warning that you try to set something you shouldn't?
    #   if things have the right size, just organize them

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
            warnings.warn(
                f"""
            No times are defined for this Rainbow.
            """
            )
        if self.nwave is None:
            warnings.warn(
                f"""
            No wavelengths are defined for this Rainbow.
            """
            )

        # does the flux have the right shape?
        if self.shape != self.flux.shape:
            message = "Flux array shape does not match (wavelength, time)."
            if self.shape == self.flux.shape[::-1]:
                warnings.warn(
                    f"""
                {message}
                Any chance it's transposed?"""
                )
            else:
                warnings.warn(
                    f"""
                {message}"""
                )

    def __repr__(self):
        """
        How should this object be represented as a string?
        """
        n = self.__class__.__name__.replace("Rainbow", "🌈")
        return f"<{n}({self.nwave}w, {self.ntime}t)>"

    def __add__(self, object):
        """
        Add the flux of a rainbow and an input array (or another rainbow)
        and output in a new rainbow object.
        Currently only supports flux addition.

        Parameters
        ----------
        object : Array or float.
            Multiple options:
            1) float
            2) 1D array with same length as wavelength axis
            3) 1D array with same length as time axis
            4) 2D array with same shape as rainbow flux
            5) Rainbow object with same dimensions as self.

        Returns
        ----------
        rainbow : Rainbow() with the same parameters as self but with added
        flux.
        """

        # Create new Rainbow() to store results in.
        result = copy.deepcopy(self)

        try:  # Object is rainbow
            if np.array_equal(
                self.wavelike["wavelength"], object.wavelike["wavelength"]
            ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
                result.fluxlike["flux"] += object.fluxlike["flux"]
            else:
                print("Objects do not share wavelength/time axes")
                return

        except AttributeError:  # Object is not rainbow

            if type(object) == float or type(object) == int:  # Float or integer.
                result.fluxlike["flux"] += object

            elif self.shape == object.shape:  # fluxlike
                result.fluxlike["flux"] += object

            elif len(object) == len(
                self.wavelike["wavelength"]
            ):  # Add along wavelength axis

                result.fluxlike["flux"] += np.transpose([object] * self.shape[1])
                if self.nwave == self.ntime:
                    raise RuntimeError(
                        f"{self} has same number of wavelengths and times; we can't tell which is which."
                    )

            elif len(object) == len(self.timelike["time"]):  # Add along time axis.
                result.fluxlike["flux"] += np.tile(object, (self.shape[0], 1))

            else:
                print("Invalid shape in addition: " + str(object.shape))
                return
        return result

    def __sub__(self, object):
        """
        Subtract the flux of a rainbow from an input array (or another rainbow)
        and output in a new rainbow object.
        Currently only supports flux subtraction.

        Parameters
        ----------
        object : Array or float.
            Multiple options:
            1) float
            2) 1D array with same length as wavelength axis
            3) 1D array with same length as time axis
            4) 2D array with same shape as rainbow flux
            5) Rainbow object with same dimensions as self.

        Returns
        ----------
        rainbow : Rainbow() with the same parameters as self but with subtracted
        flux.
        """

        # Create new Rainbow() to store results in.
        result = copy.deepcopy(self)

        try:  # Object is rainbow.
            if np.array_equal(
                self.wavelike["wavelength"], object.wavelike["wavelength"]
            ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
                result.fluxlike["flux"] -= object.fluxlike["flux"]
            else:
                print("Objects do not share wavelength/time axes")
                return

        except AttributeError:

            if type(object) == float or type(object) == int:  # Float or integer.
                result.fluxlike["flux"] -= object

            elif self.shape == object.shape:  # fluxlike
                result.fluxlike["flux"] -= object

            elif len(object) == len(
                self.wavelike["wavelength"]
            ):  # Add along wavelength axis
                result.fluxlike["flux"] -= np.transpose([object] * self.shape[1])
                if self.nwave == self.ntime:
                    raise RuntimeError(
                        f"{self} has same number of wavelengths and times; we can't tell which is which."
                    )

            elif len(object) == len(self.timelike["time"]):  # Add along time axis.
                result.fluxlike["flux"] -= np.tile(object, (self.shape[0], 1))

            else:
                print("Invalid shape in subtraction: " + str(object.shape))
                return
        return result

    def __mul__(self, object):
        """
        Multiply the flux of a rainbow and an input array (or another rainbow)
        and output in a new rainbow object.
        Currently only supports flux multiplication.

        Parameters
        ----------
        object : Array or float.
            Multiple options:
            1) float
            2) 1D array with same length as wavelength axis
            3) 1D array with same length as time axis
            4) 2D array with same shape as rainbow flux
            5) Rainbow object with same dimensions as self.

        Returns
        ----------
        rainbow : Rainbow() with the same parameters as self but with multiplied
        flux.
        """

        # Create new Rainbow() to store results in.
        result = copy.deepcopy(self)

        try:  # Object is rainbow
            if np.array_equal(
                self.wavelike["wavelength"], object.wavelike["wavelength"]
            ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
                result.fluxlike["flux"] *= object.fluxlike["flux"]
            else:
                print("Objects do not share wavelength/time axes")
                return

        except AttributeError:

            if type(object) == float or type(object) == int:  # Float or integer.
                result.fluxlike["flux"] *= object

            elif self.shape == object.shape:  # fluxlike
                result.fluxlike["flux"] *= object

            elif len(object) == len(
                self.wavelike["wavelength"]
            ):  # Add along wavelength axis
                result.fluxlike["flux"] *= np.transpose([object] * self.shape[1])
                if self.nwave == self.ntime:
                    raise RuntimeError(
                        f"{self} has same number of wavelengths and times; we can't tell which is which."
                    )

            elif len(object) == len(self.timelike["time"]):  # Add along time axis.
                result.fluxlike["flux"] *= np.tile(object, (self.shape[0], 1))

            else:
                print("Invalid shape in multiplication: " + str(object.shape))
                return
        return result

    def __truediv__(self, object):
        """
        Divide the flux of a rainbow and an input array (or another rainbow)
        and output in a new rainbow object.
        Currently only supports flux division.

        Parameters
        ----------
        object : Array or float.
            Multiple options:
            1) float
            2) 1D array with same length as wavelength axis
            3) 1D array with same length as time axis
            4) 2D array with same shape as rainbow flux
            5) Rainbow object with same dimensions as self.

        Returns
        ----------
        rainbow : Rainbow() with the same parameters as self but with divided
        flux.
        """

        # Create new Rainbow() to store results in.
        result = copy.deepcopy(self)

        try:  # Object is another rainbow.
            if np.array_equal(
                self.wavelike["wavelength"], object.wavelike["wavelength"]
            ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
                result.fluxlike["flux"] /= object.fluxlike["flux"]
            else:
                print("Objects do not share wavelength/time axes")
                return

        except AttributeError:

            if type(object) == float or type(object) == int:  # float
                result.fluxlike["flux"] /= object

            elif self.shape == object.shape:  # fluxlike array
                result.fluxlike["flux"] /= object

            elif len(object) == len(
                self.wavelike["wavelength"]
            ):  # Add along wavelength axis
                result.fluxlike["flux"] /= np.transpose([object] * self.shape[1])
                if self.nwave == self.ntime:
                    raise RuntimeError(
                        f"{self} has same number of wavelengths and times; we can't tell which is which."
                    )
            elif len(object) == len(self.timelike["time"]):  # Add along time axis.
                result.fluxlike["flux"] /= np.tile(object, (self.shape[0], 1))

            else:
                print("Invalid shape in division: " + str(object.shape))
                return
        return result

    def normalize(self):
        """
        Normalize by dividing through by the median spectrum.
        """

        # TODO, think about more careful treatment of uncertainties + good/bad data
        return self / np.nanmedian(self.flux, axis=self.timeaxis)

    def bin(self, dt=None, time=None, R=None, dw=None, wavelength=None):
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
        """

        # self.speak(f'binning')

        # bin first in time
        binned_in_time = self.bin_in_time(dt=dt, time=time)

        # then bin in wavelength
        binned = binned_in_time.bin_in_wavelength(R=R, dw=dw, wavelength=wavelength)

        # return the binned object
        return binned

    def bin_in_time(self, dt=None, time=None):
        """
        Bin the rainbow in time.

        Parameters
        ----------
        dt : astropy.units.Quantity
            The d(time) bin size for creating a grid
            that is uniform in linear space.
        time : array of astropy.units.Quantity
            An array of times, if you just want to give
            it an entirely custom array.

        The time-setting order of precendence is:
            1) time
            2) dt
        """

        # if no bin information is provided, don't bin
        if (time is None) and (dt is None):
            return self

        # set up binning parameters
        binkw = dict(weighting="inversevariance", drop_nans=False)
        # TODO: make sure drop_nans doesn't skip rows or columns unexpectedly
        if time is not None:
            binkw["newx"] = time
            # self.speak(f'binning to time={time}')
        elif dt is not None:
            binkw["dx"] = dt
            # self.speak(f'binning to dt={dt}')

        # create a new, empty Rainbow
        new = Rainbow()

        # populate the wavelength information
        new.wavelike = {**self.wavelike}
        new.metadata["wscale"] = self.wscale

        # bin the time-like variables
        # TODO (add more careful treatment of uncertainty + DQ)
        new.timelike = {}
        for k in self.timelike:
            bt, bv = bintogrid(x=self.time, y=self.timelike[k], unc=None, **binkw)
            new.timelike[k] = bv
        new.timelike[k] = bt

        # bin the flux-like variables
        # TODO (add more careful treatment of uncertainty + DQ)
        # TODO (think about cleverer bintogrid for 2D arrays)
        new.fluxlike = {}
        for k in self.fluxlike:
            # self.speak(f" binning '{k}' in time")
            # self.speak(f"  original shape was {np.shape(self.fluxlike[k])}")

            if k == "uncertainty":
                warnings.warn(
                    """
                Uncertainties aren't being handled well yet...
                """
                )
                continue

            for w in range(new.nwave):
                if self.uncertainty is None:
                    bt, bv = bintogrid(
                        x=self.time, y=self.fluxlike[k][w, :], unc=None, **binkw
                    )
                    bu = None
                else:
                    bt, bv, bu = bintogrid(
                        x=self.time,
                        y=self.fluxlike[k][w, :],
                        unc=self.uncertainty[w, :],
                        **binkw,
                    )

                if k not in new.fluxlike:
                    new_shape = (new.nwave, new.ntime)
                    new.fluxlike[k] = np.zeros(new_shape)
                    # TODO make this more robust to units

                if k == "uncertainty":
                    new.fluxlike[k][w, :] = bu
                else:
                    new.fluxlike[k][w, :] = bv
            # self.speak(f"  new shape is {np.shape(new.fluxlike[k])}")
        new._validate_core_dictionaries()
        return new

    def bin_in_wavelength(self, R=None, dw=None, wavelength=None):
        """
        Bin the rainbow in wavelength.

        Parameters
        ----------
        R : float
            The spectral resolution for creating a grid
            that is uniform in logarithmic space.
        dw : astropy.units.Quantity
            The d(wavelength) bin size for creating a grid
            that is uniform in linear space.
        wavelength : array of astropy.units.Quantity
            An array of wavelengths, if you just want to give
            it an entirely custom array.

        The wavelength-setting order of precendence is:
            1) wavelength
            2) dw
            3) R
        """

        # if no bin information is provided, don't bin
        if (wavelength is None) and (dw is None) and (R is None):
            return self

        # set up binning parameters
        binkw = dict(weighting="inversevariance", drop_nans=False)
        if wavelength is not None:
            binning_function = bintogrid
            binkw["newx"] = wavelength
            wscale = "?"
            # self.speak(f'binning to wavelength={wavelength}')
        elif dw is not None:
            binning_function = bintogrid
            binkw["dx"] = dw
            wscale = "linear"
            # self.speak(f'binning to dw={dw}')
        elif R is not None:
            binning_function = bintoR
            binkw["R"] = R
            wscale = "log"
            # self.speak(f'binning to R={R}')

        # create a new, empty Rainbow
        new = Rainbow()

        # populate the time information
        new.timelike = {**self.timelike}

        # bin the time-like variables
        # TODO (add more careful treatment of uncertainty + DQ)
        new.wavelike = {}
        for k in self.wavelike:
            bt, bv = binning_function(
                x=self.wavelength, y=self.wavelike[k], unc=None, **binkw
            )
            new.wavelike[k] = bv
        new.wavelike[k] = bt

        # bin the flux-like variables
        # TODO (add more careful treatment of uncertainty + DQ)
        # TODO (think about cleverer bintogrid for 2D arrays)
        new.fluxlike = {}
        for k in self.fluxlike:
            # self.speak(f" binning '{k}' in wavelength")
            # self.speak(f"  original shape was {np.shape(self.fluxlike[k])}")
            if k == "uncertainty":
                warnings.warn(
                    """
                Uncertainties aren't being handled well yet...
                """
                )
                continue

            for t in range(new.ntime):
                if self.uncertainty is None:
                    bt, bv = binning_function(
                        x=self.wavelength, y=self.fluxlike[k][:, t], unc=None, **binkw
                    )
                    bu = None
                else:
                    bt, bv, bu = binning_function(
                        x=self.wavelength,
                        y=self.fluxlike[k][:, t],
                        unc=self.uncertainty[:, t],
                        **binkw,
                    )

                if k not in new.fluxlike:
                    new_shape = (new.nwave, new.ntime)
                    new.fluxlike[k] = np.zeros(new_shape)
                    # TODO make this more robust to units

                if k == "uncertainty":
                    new.fluxlike[k][:, t] = bu
                else:
                    new.fluxlike[k][:, t] = bv
            # self.speak(f"  new shape is {np.shape(new.fluxlike[k])}")
        new.metadata["wscale"] = wscale
        return new

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
        Plot

        Returns
        -------
        None.

        """
        from chromatic import viz

        viz.wavelength_plot(
            self.flux,
            self.time,
            self.wavelength,
            step_size=spacing,
            ax=ax,
            plotkw={},
            fontkw={},
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
            ax.set_xlabel(
                f"Wavelength ({self.wavelength.unit.to_string('latex_inline')})"
            )
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
