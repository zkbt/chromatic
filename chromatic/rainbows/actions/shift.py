from ...imports import *

__all__ = ['shift']

def shift(self, velocity=5*u.km/u.s):
    """
    Apply a doppler shift to the wavelength array

    Parameters
    ----------
    velocity
        the systemic velocity in question, in km/s
    """
    
    new = self._create_copy()
    
<<<<<<< Updated upstream
    lightspeed = c.c.to('km/s') #speed of light in km/s
=======
    lightspeed = con.c.to('km/s') #speed of light in km/s
>>>>>>> Stashed changes
    new_wavelength = new.wavelength * ( 1 / (1 - velocity/lightspeed) )
    new.wavelike['wavelength'] = new_wavelength
    
    return new
