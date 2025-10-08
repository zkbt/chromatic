import numpy as np


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
    return flux


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
        needed to speed up computation (= batman_model
        and batman_parameters for `batman`)
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
    return flux


def exoplanet_transit(
    t,
    rp=0.1,
    per=3.0,
    t0=0.0,
    a=10.0,
    inc=90.0,
    ecc=0.0,
    w=0.0,
    u=[0.2, 0.2],
    *args,
    **kwargs,
):
    """
    One dimensional transit model with quadratic
    limb-darkening, using the `exoplanet-core` package.

    Parameters
    ----------
    t : array
        The times at which to evaluate the model.
    rp : float, array
        The radius of the planet in stellar radii
    per : float
        The period of the planet, in days.
    t0 : float
        Mid-transit time of the transit, in days.
    a : float, array
        The semi-major axis of the orbit, in stellar radii.
    inc : float, array
        The inclination of the orbit, in degrees.
    ecc: float, array
        The eccentricity of the orbit, unitless.
    w : float, array
        The argument of periastron, in degrees.
    u : array
        The quadratic limb-darkening coefficients.

    Returns
    -------
    monochromatic_flux : array
        The flux evaluated at each time.
    """

    try:
        from exoplanet_core import kepler, quad_limbdark_light_curve
    except ImportError:
        warnings.warn(
            f"""
        You're trying to produce a transit model using `exoplanet_core`,
        but it looks like you don't have `exoplanet_core` installed
        in the environment from which you're running `chromatic`.

        Please either install it with `pip install exoplanet-core`
        (see https://github.com/exoplanet-dev/exoplanet-core for details)
        or use the `.inject_transit(..., method='trapezoid')`
        option instead.
        """
        )

    # these two handy functions were adapted from the Exoplanet code:
    def warp_times(times, t_0, _pad=True):
        if _pad:
            return np.pad(t, (0, len(times)), "constant", constant_values=1) - t_0
        return times - t_0

    def get_true_anomaly(times, t_0, t_ref, n_per, e, _pad=True):
        M = (warp_times(times, t_0, _pad=_pad) - t_ref) * n_per
        if e == 0:
            return np.sin(M), np.cos(M)
        sin_f, cos_f = kepler(M, e + np.zeros(len(M)))
        return sin_f, cos_f

    n = 2 * np.pi / per
    opsw = 1 + np.sin(w)
    E0 = 2 * np.arctan2(
        np.sqrt(1 - ecc) * np.cos(w),
        np.sqrt(1 + ecc) * opsw,
    )
    M0 = E0 - ecc * np.sin(E0)
    t_periastron = t0 - M0 / n
    tref = t_periastron - t0

    # calculate the true anomaly as a function of time:
    sinf, cosf = get_true_anomaly(t, t0, tref, n, ecc, _pad=False)
    cosi = np.cos(inc * np.pi / 180)  # convert inclination to radians
    sini = np.sin(inc * np.pi / 180)  # convert inclination to radians
    w_rad = w * np.pi / 180  # convert omega to radians
    cos_w_plus_f = (np.cos(w_rad) * cosf) - (np.sin(w_rad) * sinf)  # cos(f+w)
    sin_w_plus_f = (np.sin(w_rad) * cosf) + (np.cos(w_rad) * sinf)  # sin(f+w)

    # x,y,z equations from Winn (2010):
    r = (a * (1 - (ecc**2))) / (1 + ecc * cosf)
    x = -r * cos_w_plus_f
    y = -r * sin_w_plus_f * cosi
    z = r * sin_w_plus_f * sini
    r_sky = np.sqrt(x**2 + y**2)

    # use exoplanet_core functions to extract light curve:
    # u_arr = np.array(u)
    flux = quad_limbdark_light_curve(u[0], u[1], r_sky, rp)
    # we only want the lightcurve where z > 0 (planet is between us and star)
    flux[z < 0] = 0
    return 1 + flux


transit_model_functions = dict(
    trapezoid=trapezoidal_transit, batman=batman_transit, exoplanet=exoplanet_transit
)
