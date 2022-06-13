from ...imports import *

__all__ = ['shift']

def shift(self, velocity=5*u.km/u.s):
    """
    Apply a relativistic longitudinal doppler effect to the wavelength array

    Parameters
    ----------
    velocity
        the systemic velocity in question, in km/s
    """
    
    new = self._create_copy()
    
    lightspeed = con.c.to('km/s') #speed of light in km/s

    beta = velocity/lightspeed
    new_wavelength = new.wavelength * np.sqrt( (1 - beta) / (1 + beta) )
    new.wavelike['wavelength'] = new_wavelength
    
    return new
