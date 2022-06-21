from ...imports import *

__all__ = ["shift"]


def shift(self, velocity=0 * u.km / u.s):
    """
    Doppler shift the wavelengths of this Rainbow.
    Positive velocities make wavelengths longer (redshift).
    Negative velocities make wavelengths shorter (bluesfhit).

    Parameters
    ----------
    velocity : astropy.units.Quantity
        the systemic velocity in question,
        with units of velocity (for example, u.km/u.s)
    """

    # create a new copy of this rainbow
    new = self._create_copy()

    # get the speed of light from astropy constants
    lightspeed = con.c.to("km/s")  # speed of light in km/s

    # calculate beta and make sure the units cancel
    beta = (velocity / lightspeed).decompose()

    # apply wavelength shift
    new_wavelength = new.wavelength * np.sqrt((1 + beta) / (1 - beta))
    new.wavelike["wavelength"] = new_wavelength

    # return the new object
    return new
