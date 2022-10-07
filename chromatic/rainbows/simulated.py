from .withmodel import *


class SimulatedRainbow(RainbowWithModel):
    """
    `SimulatedRainbow` objects are created from scratch
    within `chromatic`, with options for various different
    wavelength grids, time grids, noise sources, and injected
    models. They can be useful for generating quick simulated
    dataset for testing analysis and visualization tools.

    This class definition inherits from `RainbowWithModel`,
    which itself inherits from `Rainbow`.
    """

    def __init__(
        self,
        tlim=[-2.5, 2.5] * u.hour,
        dt=2 * u.minute,
        time=None,
        wlim=[0.5, 5] * u.micron,
        R=100,
        dw=None,
        wavelength=None,
        star_flux=None,
        name=None,
        signal_to_noise=None,
    ):
        """
        Initialize a `SimulatedRainbow` object from some parameters.

        This sets up an effectively empty `Rainbow` with defined
        wavelengths and times. For making more interesting
        simulated datasets, this will often be paired with
        some combination of the `.inject...` actions that inject
        various astrophysical, instrumental, or noise signatures
        into the dataset.

        The time-setting order of precendence is:
            1) time
            2) tlim + dt

        The wavelength-setting order of precendence is:
            1) wavelength
            2) wlim + dw
            3) wlim + R

        Parameters
        ----------
        tlim : list or Quantity
            The pip install -e '.[develop]'[min, max] times for creating the time grid.
            These should have astropy units of time.
        dt : Quantity
            The d(time) bin size for creating a grid
            that is uniform in linear space.
        time : Quantity
            An array of times, if you just want to give
            it an entirely custom array.
        wlim : list or Quantity
            The [min, max] wavelengths for creating the grid.
            These should have astropy units of wavelength.
        R : float
            The spectral resolution for creating a grid
            that is uniform in logarithmic space.
        dw : Quantity
            The d(wavelength) bin size for creating a grid
            that is uniform in linear space.
        wavelength : Quantity
            An array of wavelengths, if you just want to give
            it an entirely custom array.
        star_flux : numpy 1D array
            An array of fluxes corresponding to the supplied wavelengths.
            If left blank, the code assumes a normalized flux of
            flux(wavelength) = 1 for all wavelengths.
        """
        Rainbow.__init__(self)

        # (remove the history entry from creating the Rainbow)
        self._remove_last_history_entry()

        # create a history entry for this action (before other variables are defined)
        h = self._create_history_entry("SimulatedRainbow", locals())

        # set up the wavelength grid
        self._setup_fake_wavelength_grid(wlim=wlim, R=R, dw=dw, wavelength=wavelength)

        # set up the time grid
        self._setup_fake_time_grid(tlim=tlim, dt=dt, time=time)

        # save the basic inputs that aren't stored elsewhere
        self.metadata["name"] = name

        # If the flux of the star is not given,
        # assume a continuum-normlized flux where fx=1 at all wavelengths.
        if star_flux is None:
            model = np.ones(self.shape)

        # If the flux vs wavelength of the star is supplied,
        # include it in the model.
        else:
            # Check to make sure the flux and wavelengths
            # have the same shape.
            if len(star_flux) == len(self.wavelike["wavelength"]):
                model = np.transpose([star_flux] * self.shape[1])
            elif len(star_flux) == 1:
                model = star_flux * np.ones(self.shape)

        # Set uncertainty.
        self.fluxlike["flux"] = model * 1
        self.fluxlike["model"] = model * 1
        self.fluxlike["uncertainty"] = np.zeros(self.shape)

        # make sure everything is defined and sorted
        self._validate_core_dictionaries()

        if signal_to_noise is not None:
            message = f"""
            You tried to specify the noise level with
            `SimulatedRainbow(signal_to_noise={signal_to_noise})`,
            but that functionality is going away soon.
            Please replace it right now with
            `SimulatedRainbow().inject_noise(signal_to_noise={signal_to_noise})`
            so that your code will continue to work.
            You're getting away with it this time,
            but it won't work for much longer!
            """
            cheerfully_suggest(message)
            new = self.inject_noise()
            for k in ["flux", "uncertainty", "model"]:
                self.fluxlike[k] = new.fluxlike[k]

        # append the history entry to the new Rainbow
        self._record_history_entry(h)

    def _setup_fake_time_grid(
        self, tlim=[-2.5 * u.hour, 2.5 * u.hour], dt=1 * u.minute, time=None
    ):
        """
        Create a fake time grid.

        Parameters
        ----------

        tlim : list or Quantity
            The [min, max] times for creating the time grid.
            These should have astropy units of time.
        dt : Quantity
            The d(time) bin size for creating a grid
            that is uniform in linear space.
        time : Quantity
            An array of times, if you just want to give
            it an entirely custom array.

        The time-setting order of precendence is:
            1) time
            2) tlim + dt
        """
        # check we're trying to do exactly one thing
        if (tlim is None) and (time is None):
            raise RuntimeError("Please specify either `tlim` or `time`.")

        if time is None:
            t_unit = tlim[0].unit
            t_unit.to("s")
            time = np.arange(tlim[0] / t_unit, tlim[1] / t_unit, dt / t_unit) * t_unit
        else:
            t_unit = time.unit

        self.timelike["time"] = u.Quantity(time).to(u.day)
        # TODO, make this match up better with astropy time

        self._guess_tscale()

    def _setup_fake_wavelength_grid(
        self, wlim=[0.5 * u.micron, 5 * u.micron], R=100, dw=None, wavelength=None
    ):
        """
        Create a fake wavelength grid.

        Parameters
        ----------

        wlim : list or Quantity
            The [min, max] wavelengths for creating the grid.
            These should have astropy units of wavelength.
        R : float
            The spectral resolution for creating a grid
            that is uniform in logarithmic space.
        dw : Quantity
            The d(wavelength) bin size for creating a grid
            that is uniform in linear space.
        wavelength : Quantity
            An array of wavelengths, if you just want to give
            it an entirely custom array.

        The wavelength-setting order of precendence is:
            1) wavelength
            2) wlim + dw
            3) wlim + R
        """

        # check we're trying to do exactly one thing
        if (wlim is None) and (wavelength is None):
            raise RuntimeError("Please specify either `wlim` or `wavelength`.")

        # create a linear or logarithmic grid
        if wavelength is None:
            # check that we're
            if (R is None) and (dw is None):
                raise RuntimeError("Please specify either `R` or `dw`.")

            w_unit = wlim[0].unit
            if dw is None:
                self.metadata["R"] = R
                # self.metadata["wscale"] = "log"

                logw_min = np.log(wlim[0] / w_unit)
                logw_max = np.log(wlim[1] / w_unit)
                logw = np.arange(logw_min, logw_max, 1 / R)
                wavelength = np.exp(logw) * w_unit

            elif dw is not None:
                self.metadata["dw"] = dw
                # self.metadata["wscale"] = "linear"
                wavelength = (
                    np.arange(wlim[0] / w_unit, wlim[1] / w_unit, self.dw / w_unit)
                    * w_unit
                )

        # or just make sure the wavelength grid has units
        elif wavelength is not None:
            w_unit = wavelength.unit

        # make sure the wavelength array has units
        self.wavelike["wavelength"] = u.Quantity(wavelength).to(u.micron)
        self._guess_wscale()
