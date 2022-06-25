from ..imports import *

__all__ = ["calculate_planck_flux"]


def calculate_planck_flux(wavelength, temperature):
    """
    Calculate the surface flux from a thermally emitted surface,
    according to Planck function.

    Parameters
    ----------
    wavelength : u.Quantity
        The wavelengths at which to calculate,
        with units of wavelength.
    temperature : u.Quantity
        The temperature of the thermal emitter,
        with units of K.

    Returns
    -------
    surface_flux : u.Quantity
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
