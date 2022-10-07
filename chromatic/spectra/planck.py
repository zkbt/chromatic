from ..imports import *

__all__ = ["calculate_planck_flux", "get_planck_photons"]


def calculate_planck_flux(wavelength, temperature):
    """
    Calculate the surface flux from a thermally emitted surface,
    according to Planck function.

    Parameters
    ----------
    wavelength : Quantity
        The wavelengths at which to calculate,
        with units of wavelength.
    temperature : Quantity
        The temperature of the thermal emitter,
        with units of K.

    Returns
    -------
    surface_flux : Quantity
        The surface flux, evaluated at the wavelengths.
    """

    # define variables as shortcut to the constants we need
    h = con.h
    k = con.k_B
    c = con.c

    # the thing that goes in the exponent (its units better cancel!)
    z = h * c / (wavelength * k * temperature)

    # calculate the intensity from the Planck function
    intensity = (2 * h * c**2 / wavelength**5 / (np.exp(z) - 1)) / u.steradian

    # calculate the flux assuming isotropic emission
    flux = np.pi * u.steradian * intensity

    # return the intensity
    return flux.to("W/(m**2*nm)")


def get_planck_photons(
    temperature=3000, wavelength=None, R=100, wlim=[0.04, 6] * u.micron, **kw
):
    """
    Calculate the surface flux from a thermally emitted surface,
    according to Planck function, in units of photons/(s * m**2 * nm).

    Parameters
    ----------
    temperature : Quantity
        The temperature of the thermal emitter,
        with units of K.
    wavelength : Quantity, optional
        The wavelengths at which to calculate,
        with units of wavelength.
    R : float, optional
        The spectroscopic resolution for creating a log-uniform
        grid that spans the limits set by `wlim`, only if
        `wavelength` is not defined.
    wlim : Quantity, optional
        The two-element [lower, upper] limits of a wavelength
        grid that would be populated with resolution `R`, only if
        `wavelength` is not defined.
    **kw : dict, optional
        Other keyword arguments will be ignored.

    Returns
    -------
    photons : Quantity
        The surface flux in photon units

    This evaluates the Planck function at the exact
    wavelength values; it doesn't do anything fancy to integrate
    over binwidths, so if you're using very wide (R~a few) bins
    your integrated fluxes will be messed up.

    """

    # make sure the temperature unit is good (whether or not it's supplied)
    temperature_unit = u.Quantity(temperature).unit
    if temperature_unit == u.K:
        temperature_with_unit = temperature
    elif temperature_unit == u.Unit(""):
        temperature_with_unit = temperature * u.K

    # create a wavelength grid if one isn't supplied
    if wavelength is None:
        wavelength_unit = wlim.unit
        wavelength = (
            np.exp(np.arange(np.log(wlim[0].value), np.log(wlim[1].value), 1 / R))
            * wavelength_unit
        )

    energy = calculate_planck_flux(
        wavelength=wavelength, temperature=temperature_with_unit
    )
    photon_energy = con.h * con.c / wavelength / u.ph

    return wavelength, (energy / photon_energy).to(u.ph / (u.s * u.m**2 * u.nm))
