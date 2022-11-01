from ...imports import *
import batman

__all__ = ["inject_transit"]


def trapezoidal_transit(t, delta=0.01, P=1, t0=0, T=0.1, tau=0.01, baseline=1.0):

    """
    One dimensional Trapezoid Transit model,
    using the symbols defined for a circular
    transit approximation in Winn (2010).

    This is a fittable astropy model.

    Parameters
    ----------
    delta : float
        The depth of the transit, as a fraction of the out-of-transit flux.
    P : float
        The period of the planet, in days.
    t0 : float
        Mid-transit time of the transit, in days.
    T : float
        The duration of the transit (from mid-ingress to mid-egress), in days.
    tau : float
        The duration of ingress/egress, in days.
    baseline : float
        The baseline, out-of-transit flux level.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from henrietta.fitting import TrapezoidTransit

        plt.figure()
        model = TrapezoidTransit()
        t = np.arange(-5, 5, .01)

        for x in np.linspace(0.02, 0.2, 4):
            model.delta = x
            model.T = x
            plt.plot(t, model(t), lw=2)

        plt.show()
    """

    # calculate a phase-folded time (still in units of days)
    x = (t - t0 + 0.5 * P) % P - 0.5 * P

    # Compute the four points where the trapezoid changes slope
    # x1 <= x2 <= x3 <= x4

    if tau > T:
        x1 = -tau
        x2 = 0
        x3 = 0
        x4 = tau
    else:
        x1 = -(T + tau) / 2.0
        x2 = -(T - tau) / 2.0
        x3 = (T - tau) / 2.0
        x4 = (T + tau) / 2.0

    # Compute model values in pieces between the change points
    range_a = np.logical_and(x >= x1, x < x2)
    range_b = np.logical_and(x >= x2, x < x3)
    range_c = np.logical_and(x >= x3, x < x4)

    if tau == 0:
        slope = np.inf
    else:
        slope = delta / tau
    val_a = 1 - slope * (x - x1)
    val_b = 1 - delta
    val_c = 1 - slope * (x4 - x)
    result = (
        np.select([range_a, range_b, range_c], [val_a, val_b, val_c], default=1)
        * baseline
    )
    return result


def exoplanetcore_transit(*args, **kwargs):
    raise NotImplementedError(
        "The `exoplanet-core` transit option isn't available yet. Sorry!"
    )


def inject_transit(
    self, depth=0.01, radius_ratio=None, method="trapezoid", **transit_parameters
):

    """
    Simulate a wavelength-dependent planetary transit.

    This uses one of a few methods to inject a transit
    signal into the `Rainbow`, allowing the transit
    depth to change with wavelength (for example due to a
    planet's effective radius changing with wavelength due
    to its atmospheric transmission spectrum). Other
    parameters can also be wavlength-dependent, but
    some (like period, inclination, etc...) probably
    shouldn't be.

    The current methods include:

    `'trapezoid'` to inject a cartoon transit, using nomenclature
    from [Winn (2010)](https://arxiv.org/abs/1001.2010).
    This is the default method, to avoid package dependencies
    that can be finicky to compile and/or install on different
    operating systems.

    `'batman'` to inject a limb-darkened transit using [batman-package](https://lkreidberg.github.io/batman/docs/html/index.html)
    This method requires that `batman-package` be installed,
    and it will try to throw a helpful warning message if
    it's not.

    `'exoplanet-core'` to inject a limb-darkened transit using [exoplanet-core](https://github.com/exoplanet-dev/exoplanet-core).
    This option is still under developed. Coming soon! We hope
    this will be the best of both worlds, prodiving a limb-darkened
    transit without scary installation dependencies.

    Parameters
    ----------
    depth : float, array, None
        The transit depth = [(planet radius)/(star radius)]**2,
        which can be either a single value for all wavelengths,
        or an array with one value for each wavelength.
        Only one of [`depth`, `radius_ratio`] can be not None.
        (default 0.01)
    radius_ratio : float, array, None
        The planet-to-star radius ratio = [transit depth]**0.5,
        which can be either a single value for all wavelengths,
        or an array with one value for each wavelength.
        Only one of [`depth`, `radius_ratio`] can be not None.
        (default None)
    method : str
        What method should be used to inject the transits? Different
        methods will produce different results and have different options.
        The currently implement options are `'trapezoid'` and `'batman'`.
    **transit_parameters : dict
        All additional keywords will be passed to the transit model.
        The accepted keywords for the different methods are as follows.
            `'trapezoid'` accepts the following keyword arguments:
                "delta" = The depth of the transit, as a fraction of the out-of-transit flux (default 0.01)
                "P" = The orbital period of the planet, in days (default 1.0)
                "t0" = Mid-transit time of the transit, in days (default 0.0)
                "T" = The duration of the transit (from mid-ingress to mid-egress), in days (default 0.1)
                "tau" = The duration of ingress/egress, in days (default 0.01)
                "baseline" = The baseline, out-of-transit flux level (default 1.0)
            `batman` accepts the following keyword arguments:
                "rp" = (planet radius)/(star radius), unitless (default 0.1)
                "t0" = Mid-transit time of the transit, in days (default 0.0)
                "per" = The orbital period of the planet, in days (default 1.0)
                "a" = (semi-major axis)/(star radius), unitless (default 10)
                "inc" = The orbital inclination, in degrees (default 90)
                "ecc" = The orbital eccentricity, unitless (default 0.0)
                "w" = The longitude of periastron, in degrees (default 0.0)
                "limb_dark" = The limb-darkening model (default "quadratic"), possible
                    values described in more detail in batman documentation.
                "u" = The limb-darkening coefficients (default [0.2, 0.2])
                    These coefficients can be:
                        -one value (if limb-darkening law requires only one value)
                        -a 1D list/array of coefficients for constant limb-darkening
                        -a 2D array of the form (n_wavelengths, n_coefficients) where
                        each row is the set of limb-darkening coefficients corresponding
                        to a single wavelength
                    Note that this currently does not calculate the appropriate
                    coefficient vs wavelength variations itself-there exist codes
                    (such as hpparvi/PyLDTk and nespinoza/limb-darkening) which
                    can be used for this.
    """

    '''# create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inject_transit", locals())
    h = h.replace("transit_parameters={", "**{")

    if ((depth == None) and (radius_ratio == None)) or ((depth ! None) and (radius_ratio != None)):
        raise ValueError(
            f"""
        You appear to have called `inject_transit` with
        the following inputs:
            depth={depth}
            radius_ratio={radius_ratio}
        Exactly one of these needs to be set,
        and the other needs to be None.
        """
        )

    # create a copy of the existing Rainbow
    new = self._create_copy()

    # Defaults for planet simulation.
    if method == "trapezoid":
        defaults = {
            "delta": 0.01,
            "P": 1.0,
            "t0": 0.0,
            "T": 0.1,
            "tau": 0.01,
            "baseline": 1.0,
        }
    elif method == "batman":
        defaults = {
            "t0": 0.0,
            "per": 1.0,
            "a": 10.0,
            "inc": 90.0,
            "ecc": 0.0,
            "w": 0.0,
            "limb_dark": "quadratic",
            "u": [0.2, 0.2],
        }

    # first, make sure depth has the right dimensions
    if type(planet_radius) != float and len(planet_radius) != self.nwave:
        print(
            "Invalid planet radius array: must be float or have shape "
            + str(np.shape(self.wavelike["wavelength"]))
        )

    # Read in planet parameters.
    for i in range(len(transit_parameters.keys())):
        key = list(transit_parameters.keys())[i]
        if key in list(defaults.keys()):
            defaults[key] = transit_parameters[key]
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
            m = batman.TransitModel(params, new.time.to_value("day"))
        planet_flux[i] = m.light_curve(params)

    new.planet_model = planet_flux
    new.flux *= new.planet_model
    new.model = new.fluxlike.get("model", 1) * new.planet_model

    parameters = dict(**defaults)
    parameters.update(rp_unbinned=rprs, u_unbinned=u_arr)
    new.metadata["transit_parameters"] = parameters
    new.metadata["transit_version"] = f"batman-package=={batman.__version__}"
    if np.ndim(u) == 1:
        new.wavelike["injected-u"] = u_arr
    else:
        for i in range(np.shape(u_arr)[1]):
            new.wavelike[f"injected-transit-u{i+1}"] = u_arr[:, i]
    new.wavelike["injected-transit-rp"] = rprs
    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new object
    return new
'''
