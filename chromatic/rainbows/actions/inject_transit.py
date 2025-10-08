from ...tools.transits import *

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


def inject_transit(
    self,
    planet_radius=0.1,
    method="exoplanet",
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

    `'exoplanet'` to inject a limb-darkened transit using [exoplanet-core](https://github.com/exoplanet-dev/exoplanet-core).
    This option requires `exoplanet-core` be installed,
    but it doesn't require complicated dependencies or
    compiling steps, so it's already included as a dependency.

    Parameters
    ----------
    planet_radius : float, array, None
        The planet-to-star radius ratio = [transit depth]**0.5,
        which can be either a single value for all wavelengths,
        or an array with one value for each wavelength.
    method : str
        What method should be used to inject the transits? Different
        methods will produce different results and have different options.
        The currently implement options are `'trapezoid'` and `'exoplanet'`.
    **transit_parameters : dict
        All additional keywords will be passed to the transit model.
        The accepted keywords for the different methods are as follows.
            `'trapezoid'` accepts the following keyword arguments:
                `delta` = The depth of the transit, as a fraction of the out-of-transit flux (default 0.01)
                (If not provided, it will be set by `planet_radius`.)
                `P` = The orbital period of the planet, in days (default 3.0)
                `t0` = Mid-transit time of the transit, in days (default 0.0)
                `T` = The duration of the transit (from mid-ingress to mid-egress), in days (default 0.1)
                `tau` = The duration of ingress/egress, in days (default 0.01)
                `baseline` = The baseline, out-of-transit flux level (default 1.0)
            `'exoplanet'` accepts the following keyword arguments:
                `rp` = (planet radius)/(star radius), unitless (default 0.1)
                (If not provided, it will be set by `planet_radius`.)
                `t0` = Mid-transit time of the transit, in days (default 0.0)
                `per` = The orbital period of the planet, in days (default 3.0)
                `a` = (semi-major axis)/(star radius), unitless (default 10)
                `inc` = The orbital inclination, in degrees (default 90)
                `ecc` = The orbital eccentricity, unitless (default 0.0)
                `w` = The longitude of periastron, in degrees (default 0.0)
                `u` = The quadratic limb-darkening coefficients (default [0.2, 0.2])
                    These coefficients can only be a 2D array of the form (n_wavelengths, n_coefficients) where
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
    elif method == "exoplanet":
        parameters_to_use = {
            "rp": planet_radius,
            "t0": 0.0,
            "per": 3.0,
            "a": 10.0,
            "inc": 90.0,
            "ecc": 0.0,
            "w": 0.0,
            "u": [[0.2, 0.2]],
        }
    else:
        raise ValueError(
            f"""
        'method' must be one of ['exoplanet', 'trapezoid']
        """
        )

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
    planet_flux = np.ones(new.shape)
    for i in range(self.nwave):
        parameters_for_this_wavelength = {
            k: get_for_wavelength(parameters_to_use[k], i) for k in parameters_to_use
        }
        f = transit_model_functions[method]
        monochromatic_flux = f(t, **parameters_for_this_wavelength)
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
