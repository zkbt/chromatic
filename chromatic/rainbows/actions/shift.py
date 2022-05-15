from ...imports import *

def trim(self, velocity=5*u.km/u.s):
    """
    Apply a doppler shift to the wavelength array

    Parameters
    ----------
    velocity
        the systemic velocity in question, in km/s
    """
    c = 3e5 * u.km / u.s #speed of light in km/s
    shifted = self.wave * ( 1 / (1 + v/c) )
    
    return shifted
