from .rainbow import *


class SimulatedRainbow(Rainbow):
    def __init__(
        self,
        signal_to_noise=100,
        tlim=[-2.5, 2.5] * u.hour,
        dt=2 * u.minute,
        time=None,
        wlim=[0.5, 5] * u.micron,
        R=100,
        dw=None,
        wavelength=None,
        star_flux=None,
        name=None,
    ):
        """
        Create a simulated rainbow object.

        Parameters
        ----------

        signal_to_noise : float
            The signal-to-noise per wavelength per time.
            For example, S/N=100 would mean that the
            uncertainty on the flux for each each
            wavelength-time data point will be 1%.

        tlim : list or array of astropy.units.Quantity
            The pip install -e '.[develop]'[min, max] times for creating the time grid.
            These should have astropy units of time.
        dt : astropy.units.Quantity
            The d(time) bin size for creating a grid
            that is uniform in linear space.
        time : array of astropy.units.Quantity
            An array of times, if you just want to give
            it an entirely custom array.

        The time-setting order of precendence is:
            1) time
            2) tlim + dt

        wlim : list or array of astropy.units.Quantity
            The [min, max] wavelengths for creating the grid.
            These should have astropy units of wavelength.
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
            2) wlim + dw
            3) wlim + R

        star_flux : numpy 1D array
            An array of fluxes corresponding to the supplied wavelengths.
            If left blank, the code assumes a normalized flux of
            flux(wavelength) = 1 for all wavelengths.
        planet : boolean
            Set to True to add a transiting planet.
        planet_params: Dictionary
            Stores planetary parameters to model, read description for
            Rainbow.add_planet transit() for more details.
        planet_radius: float or 1D numpy array
            Planet Radius/Stellar Radius of planet.  Use 1D if wavelength-
            dependent.  Must be same shape as wavelength array.
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

        # Save SNR.
        self.metadata["signal_to_noise"] = signal_to_noise

        # If the flux of the star is not given,
        # assume a continuum-normlized flux where fx=1 at all wavelengths.
        if star_flux is None:
            self.fluxlike["model"] = np.ones(self.shape)

        # If the flux vs wavelength of the star is supplied,
        # include it in the model.
        else:
            # Check to make sure the flux and wavelengths
            # have the same shape.
            if len(star_flux) == len(self.wavelike["wavelength"]):
                self.fluxlike["model"] = np.transpose([star_flux] * self.shape[1])

        # Set uncertainty.
        self.fluxlike["uncertainty"] = self.fluxlike["model"] / signal_to_noise
        self.fluxlike["flux"] = np.random.normal(
            self.fluxlike["model"], self.fluxlike["uncertainty"]
        )

        # append the history entry to the new Rainbow
        self._record_history_entry(h)
        self.metadata["name"] = name

    def _setup_fake_time_grid(
        self, tlim=[-2.5 * u.hour, 2.5 * u.hour], dt=1 * u.minute, time=None
    ):
        """
        Create a fake time grid.

        Parameters
        ----------

        tlim : list or array of astropy.units.Quantity
            The [min, max] times for creating the time grid.
            These should have astropy units of time.
        dt : astropy.units.Quantity
            The d(time) bin size for creating a grid
            that is uniform in linear space.
        time : array of astropy.units.Quantity
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

    def _setup_fake_wavelength_grid(
        self, wlim=[0.5 * u.micron, 5 * u.micron], R=100, dw=None, wavelength=None
    ):
        """
        Create a fake wavelength grid.

        Parameters
        ----------

        wlim : list or array of astropy.units.Quantity
            The [min, max] wavelengths for creating the grid.
            These should have astropy units of wavelength.
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
