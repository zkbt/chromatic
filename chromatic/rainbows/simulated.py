from .rainbow import *
import batman as batman


class SimulatedRainbow(Rainbow):
    def __init__(
        self,
        signal_to_noise=100,
        tlim=[-2.5, 2.5] * u.hour,
        dt=10 * u.minute,
        time=None,
        wlim=[0.5, 5] * u.micron,
        R=20,
        dw=None,
        wavelength=None,
        star_flux=None,
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

        self.timelike["time"] = u.Quantity(time)
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
        self.wavelike["wavelength"] = u.Quantity(wavelength)
        self._guess_wscale()

        # this should break if the units aren't length
        w_unit.to("m")

    def inject_transit(self, planet_params={}, planet_radius=0.1):

        """
        Simulate a wavelength-dependent planetary transit using
        batman.

        Parameters
        ----------

        planet_params : Dictionary
            Values for planetary parameters for use in batman modelling.
            Any values not supplied will be set to defaults:
                "t0" = time of inferior conjunction (days) (default 0)
                "per" = orbital period (days) (detault 1)
                "a" = semi-major axis (units of stellar radii) (default 15)
                "inc" = inclination (degrees) (default 90)
                "ecc" = eccentricity (default 0)
                "w" = longitude of periastron (degrees)(default 0)
                "limb_dark" = limb-darkening model (default "nonlinear"), possible
                    values described in more detail in batman documentation
                "u" = limb-darkening coefficients (default [0.5, 0.1, 0.1, -0.1])
                    Can take 3 forms:
                    -A single value (if limb-darkening law requires only one value)
                    -A 1D list/array of coefficients corresponding to the limb-darkening
                    law
                    -A 2D array of the form (n_wavelengths, n_coefficients) where
                    each row is the set of limb-darkening coefficients corresponding
                    to a single wavelength

                    Note that this currently does not calculate the appropriate
                    coefficient vs wavelength variations itself- there exist codes
                    (such as hpparvi/PyLDTk and nespinoza/limb-darkening) which
                    can be used for this.

            example value: planet_params = {"a":12, "inc":87}

        planet_radius = Two options:
                1D array with same dimensions as wavelength array,
                    each value corresponds to planet radius/stellar radius at that
                    wavelength.
                float representing Rp/Rstar if the radius is not wavelength-dependent.
            example value: planet_radius = 0.01,

        """

        # First, make sure planet_radius has the right dimension.
        if type(planet_radius) != float and len(planet_radius) != self.nwave:
            print(
                "Invalid planet radius array: must be float or have shape "
                + str(np.shape(self.wavelike["wavelength"]))
            )

        # Defaults for planet simulation.
        defaults = {
            "t0": 0,
            "per": 3,
            "a": 10,
            "inc": 90,
            "ecc": 0,
            "w": 0,
            "limb_dark": "nonlinear",
            "u": [0.5, 0.1, 0.1, -0.1],
        }

        # Read in planet parameters.
        for i in range(len(planet_params.keys())):
            key = list(planet_params.keys())[i]
            if key in list(defaults.keys()):
                defaults[key] = planet_params[key]
            else:
                print("Warning: " + str(key) + " not a valid parameter")

        # Initialize batman model.
        params = batman.TransitParams()
        params.t0 = defaults["t0"]
        params.per = defaults["per"]
        params.a = defaults["a"]
        params.inc = defaults["inc"]
        params.ecc = defaults["ecc"]
        params.w = defaults["w"]
        params.limb_dark = defaults["limb_dark"]

        # Deal with limb-darkening.
        if len(np.shape(defaults["u"])) < 2:  # Coefficients constant with wavelength
            u_arr = np.tile(defaults["u"], (self.nwave, 1))

        elif (
            len(np.shape(defaults["u"])) == 2
        ):  # 2D array of coefficients, along wavelength axis
            if np.shape(defaults["u"])[0] != self.nwave:
                print("Shape of limb-darkening array does not match wavelengths.")
                return
            u_arr = defaults["u"]

        else:
            print("Invalid limb-darkening coefficient array.")

        # Read in planetary radius.
        if type(planet_radius) == float:
            rprs = np.zeros(self.nwave) + planet_radius
        else:
            rprs = planet_radius

        planet_flux = np.zeros((self.nwave, self.ntime))

        for i in range(self.nwave):
            params.rp = rprs[i]
            params.u = u_arr[i]
            # print(params.u)
            try:
                m
            except NameError:
                m = batman.TransitModel(params, self.timelike["time"].to("day").value)
            planet_flux[i] = m.light_curve(params)

        result = self._create_copy()
        result.fluxlike["model"] *= planet_flux
        result.fluxlike["flux"] *= planet_flux

        return result
