from ...imports import *

def shift(self, velocity=5*u.km/u.s):
    """
    Apply a doppler shift to the wavelength array

    Parameters
    ----------
    velocity
        the systemic velocity in question, in km/s
    """
    
    new = self._create_copy()
    
    c = 3e5 * u.km / u.s #speed of light in km/s
    new_wavelength = new.wavelength * ( 1 / (1 + v/c) )
    new.wavelike['wavelength'] = new_wavelength
    
    return new
