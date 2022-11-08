from ...imports import *

__all__ = ["inject_transit"]


def get_for_wavelength(x, i=0):
    """
    Get a parameter value, which might be wavelength-specific.

    This is a simple wrapper to help handle wavelength-dependent
    and wavelength-shared parameters in the same way.

    Parameters
    ----------
    x : float, array, str
        Some parameter going into a model.
    i : int
        An index, corresponding to a wavelength.

    Returns
    -------
    some_parameter : float, array, str
        The parameter value(s) for a particular wavelength,
        which might have been shared across all wavelengths
        or might be unique to just this wavelength.
    """
    # assume strings are constant across wavelength
    if type(x) == str:
        return x
    # don't try to index scalars
    if np.shape(x) == ():
        return x
    # assume (1, ?) arrays are constant across wavelength
    if np.shape(x)[0] == 1:
        return x[0]
    # extract for a particular wavelength
    if np.shape(x)[0] > 1:
        return x[i]


def trapezoidal_transit(t, delta=0.01, P=1, t0=0, T=0.1, tau=0.01, baseline=1.0):
    """
    One dimensional Trapezoid Transit model,
    using the symbols defined for a circular
    transit approximation in Winn (2010).

    Parameters
    ----------
    t : array
        The times at which to evaluate the model.
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

    Returns
    -------
    monochromatic_flux : array
        The flux evaluated at each time.
    cached_inputs : dict
        A kludge to store any intermediate variables
        needed to speed up computation (= none for
        the trapezoidal transit).
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
    flux = (
        np.select([range_a, range_b, range_c], [val_a, val_b, val_c], default=1)
        * baseline
    )
    return flux, {}


def batman_transit(
    t,
    rp=0.1,
    t0=0.0,
    per=3.0,
    a=10.0,
    inc=90.0,
    ecc=0.0,
    w=0.0,
    u=[0.2, 0.2],
    limb_dark="quadratic",
    batman_model=None,
    batman_params=None,
    **kw,
):
    """
    One dimensional batman Transit model,
    using variables defined with `batman-package`.

    Parameters
    ----------
    t : array
        The times at which to evaluate the model.
    rp : float, array
        The radius of the planet, in stellar radii.
    t0 : float, array
        Mid-transit time of the transit, in days.
    per : float, array
        The period of the planet, in days.
    a : float, array
        The semi-major axis of the orbit, in stellar radii.
    inc : float, array
        The inclination of the orbit, in degrees.
    ecc: float, array
        The eccentricity of the orbit, unitless.
    w : float, array
        The argument of periastron, in degrees.
    limb_dark : str
        The limb-darkening law to use.
    u : array
        The limb-darkening coefficients, in a 2D array of
        shape (nwavelengths, ncoefficients_per_wavelength).
        For example, the coefficients for wavelength-shared
        quadratic limb-darkening might have a shape of (1, 2).
    **kw : dict
        All additional parameters will be passed to the
        batman.TransitModel initialization. Possibilities
        with defaults include [max_err=1.0, nthreads=1, fac=None,
        transittype='primary', supersample_factor=1, exp_time=0.0]


    Returns
    -------
    monochromatic_flux : array
        The flux evaluated at each time.
    cached_inputs : dict
        A kludge to store any intermediate variables
        needed to speed up computation (= none for
        the trapezoidal transit).
    """

    try:
        import batman
    except ImportError:
        warnings.warn(
            f"""
        You're trying to produce a transit model using `batman`,
        but it looks like you don't have `batman-package` installed
        in the environment from which you're running `chromatic`.

        Please either install it with `pip install batman-package`
        (see https://github.com/lkreidberg/batman for details)
        or use the `.inject_transit(..., method='trapezoid')`
        option instead.
        """
        )

    # set the batman parameters
    if batman_params is None:
        batman_params = batman.TransitParams()
    batman_params.rp = rp
    batman_params.t0 = t0
    batman_params.per = per
    batman_params.a = a
    batman_params.inc = inc
    batman_params.ecc = ecc
    batman_params.w = w
    batman_params.u = u
    batman_params.limb_dark = limb_dark

    # generate the batman model
    if batman_model is None:
        batman_model = batman.TransitModel(batman_params, t)

    flux = batman_model.light_curve(batman_params)
    cached_inputs = dict(batman_model=batman_model, batman_params=batman_params)
    return flux, cached_inputs


def exoplanet_transit(*args, **kwargs):
    raise NotImplementedError(
        "The `exoplanet-core` transit option isn't available yet. Sorry!"
    )


transit_model_functions = dict(
    trapezoid=trapezoidal_transit, batman=batman_transit, exoplanet=exoplanet_transit
)


def inject_transit(
    self,
    planet_radius=0.1,
    method="trapezoid",
    **transit_parameters,
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
    planet_radius : float, array, None
        The planet-to-star radius ratio = [transit depth]**0.5,
        which can be either a single value for all wavelengths,
        or an array with one value for each wavelength.
    method : str
        What method should be used to inject the transits? Different
        methods will produce different results and have different options.
        The currently implement options are `'trapezoid'` and `'batman'`.
    **transit_parameters : dict
        All additional keywords will be passed to the transit model.
        The accepted keywords for the different methods are as follows.
            `'trapezoid'` accepts the following keyword arguments:
                `delta` = The depth of the transit, as a fraction of the out-of-transit flux (default 0.01)
                (If not provided, it will be set by `depth` or `radius_ratio`.)
                `P` = The orbital period of the planet, in days (default 3.0)
                `t0` = Mid-transit time of the transit, in days (default 0.0)
                `T` = The duration of the transit (from mid-ingress to mid-egress), in days (default 0.1)
                `tau` = The duration of ingress/egress, in days (default 0.01)
                `baseline` = The baseline, out-of-transit flux level (default 1.0)
            `'batman'` accepts the following keyword arguments:
                `rp` = (planet radius)/(star radius), unitless (default 0.1)
                (If not provided, it will be set by depth or radius_ratio.)
                `t0` = Mid-transit time of the transit, in days (default 0.0)
                `per` = The orbital period of the planet, in days (default 1.0)
                `a` = (semi-major axis)/(star radius), unitless (default 10)
                `inc` = The orbital inclination, in degrees (default 90)
                `ecc` = The orbital eccentricity, unitless (default 0.0)
                `w` = The longitude of periastron, in degrees (default 0.0)
                `limb_dark` = The limb-darkening model (default "quadratic"), possible
                    values described in more detail in batman documentation.
                `u` = The limb-darkening coefficients (default [0.2, 0.2])
                    These coefficients can be:
                        -one value (if limb-darkening law requires only one value)
                        -a 1D list/array of coefficients for constant limb-darkening
                        -a 2D array of the form (n_wavelengths, n_coefficients) where
                        each row is the set of limb-darkening coefficients corresponding
                        to a single wavelength
                    Note that this currently does not calculate the appropriate
                    coefficient vs wavelength variations itself; there exist codes
                    (such as hpparvi/PyLDTk and nespinoza/limb-darkening) which
                    can be used for this.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inject_transit", locals())
    h = h.replace("transit_parameters={", "**{")

    # create a copy of the existing Rainbow
    new = self._create_copy()

    # make sure the depth is set, with some flexibility
    # to allow for different names. parameter names that
    # belong directly to the transit model [delta, rp]
    # will take precendence first, then [depth], then
    # [planet_radius = the default]

    # set defaults for planet simulation
    if method == "trapezoid":
        parameters_to_use = {
            "delta": planet_radius**2 * np.sign(planet_radius),
            "P": 1.0,
            "t0": 0.0,
            "T": 0.1,
            "tau": 0.01,
            "baseline": 1.0,
        }
    elif method == "batman":
        parameters_to_use = {
            "rp": planet_radius,
            "t0": 0.0,
            "per": 3.0,
            "a": 10.0,
            "inc": 90.0,
            "ecc": 0.0,
            "w": 0.0,
            "limb_dark": "quadratic",
            "u": [[0.2, 0.2]],
        }

    # update based on explicit keyword arguments
    parameters_to_use.update(**transit_parameters)

    # check the parameter shapes are legitimate
    for k, v in parameters_to_use.items():
        s = np.shape(v)
        if (s != ()) and (s[0] not in [1, new.nwave]):
            raise ValueError(
                f"""
            The parameter {k}={v}
            has a shape of {np.shape(v)}, which we don't know
            how to interpret. It should be a single value,
            or have a first dimension of either 1 or nwave={new.nwave}.
            """
            )

    # call the model for each wavelength
    t = new.time.to_value("day")
    cached_inputs = {}
    planet_flux = np.ones(new.shape)
    for i in range(self.nwave):
        parameters_for_this_wavelength = {
            k: get_for_wavelength(parameters_to_use[k], i) for k in parameters_to_use
        }
        f = transit_model_functions[method]
        monochromatic_flux, cached_inputs = f(
            t, **parameters_for_this_wavelength, **cached_inputs
        )
        planet_flux[i, :] = monochromatic_flux

    # store the model in the new Rainbow object
    new.planet_model = planet_flux
    new.flux *= new.planet_model
    new.model = new.fluxlike.get("model", 1) * new.planet_model

    # store the injected parameters as metadata or wavelike
    new.metadata["injected_transit_method"] = method
    new.metadata["injected_transit_parameters"] = parameters_to_use
    for k, v in parameters_to_use.items():
        label = f"injected_transit_{k}"
        s = np.shape(v)
        if s == ():
            continue
        elif s[0] == new.nwave:
            if len(s) == 1:
                new.wavelike[label] = v
            elif len(s) > 1:
                for i in range(s[1]):
                    new.wavelike[f"{label}{i+1}"] = v[:, i]

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new object
    return new
