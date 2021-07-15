from ..imports import *
from ..talker import Talker
from ..resampling import *
from .readers import *

class Rainbow(Talker):
    """
    Rainbow objects represent the flux of an object as a function of both wavelength and time.
    """

    _core_dictionaries = ["fluxlike", "timelike", "wavelike", "metadata"]

    # define which axis is which
    waveaxis = 0
    timeaxis = 1

    def __init__(self, wavelength=None, time=None, flux=None, uncertainty=None, **kw):
        """
        Initialize a Rainbow object.

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

        # metadata are arbitrary types of information we need
        self.metadata = {}

        # wavelike quanities are 1D arrays with nwave elements
        self.wavelike = {}
        self.wavelike["wavelength"] = wavelength

        # timelike quantities are 1D arrays with ntime elements
        self.timelike = {}
        self.timelike["time"] = time

        # fluxlike quantities are 2D arrays with nwave x time elements
        self.fluxlike = {}
        self.fluxlike["flux"] = flux
        self.fluxlike["uncertainty"] = uncertainty

        self._guess_wscale()

    def _guess_wscale(self):
        """
        Try to guess the wscale from the wavelengths.
        """

        # give up if there's no wavelength array
        if self.wavelength is None:
            return "?"

        # calculate difference arrays
        w = self.wavelength.value
        dw = np.diff(w)
        dlogw = np.diff(np.log(w))

        # test the three options
        if np.all(dw == dw[0]):
            self.metadata["wscale"] = "linear"
        elif np.all(dlogw == dlogw[0]):
            self.metadata["wscale"] = "log"
        else:
            self.metadata["wscale"] = "?"

    @property
    def wavelength(self):
        return self.wavelike['wavelength']

    @property
    def time(self):
        return self.timelike['time']

    @property
    def flux(self):
        return self.fluxlike['flux']

    @property
    def uncertainty(self):

        return self.fluxlike.get('uncertainty', None)

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
        The shape of the flux array (nwave, ntime)
        """
        return (self.nwave, self.ntime)

    @property
    def nwave(self):
        """
        The number of wavelengths
        """
        if self.wavelength is None:
            return 0
        else:
            return len(self.wavelength)

    @property
    def ntime(self):
        """
        The number of times
        """
        if self.time is None:
            return 0
        else:
            return len(self.time)

    @property
    def nflux(self):
        """
        The number of fluxes
        """
        return np.prod(self.shape)

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
        '''
        Normalize by dividing through by the median spectrum.
        '''

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

            if k == 'uncertainty':
                warnings.warn("Uncertainties aren't being handled well yet...")
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
            if k == 'uncertainty':
                warnings.warn("Uncertainties aren't being handled well yet...")
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
            ax = plt.gca()

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

Rainbow.from_eureka = from_eureka
Rainbow.from_x1dints = from_x1dints
