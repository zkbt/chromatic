from ...imports import *

__all__ = ['shift']

from chromatic import *

# create an example object to test on 
tester = SimulatedRainbow()

# define your new method - calculate delta D_spot (lambda)
def starspot_depth(self,
                          T_spot = 4000 *u.K,
                          T_unspot = 3500*u.K,
                          f_spot = 0.1,
                          f_tra = 0.,
                          r_sun = 0.5 * u.Rsun,
                          m_sun = 0.5 * u.Msun, ):
    """
    Calculate delta D(lambda), the contribution of 
    starspot features to the wavelength dependent transit depth

    Parameters
    ----------
    T1
        the photosphere temperature, in K
    T2
        the spot temperature, in K
    f_spot
        the global spot coverage fraction
    f_tra
        the transit chord spot coverage fraction
    r_sun
        the stellar radius, in solar radii
    m_sun
        the stellar mass, in solar masses
    """
    
    new = self._create_copy()
    
    S_spot = Star(teff=T_spot, radius=r_sun,mass=m_sun,R=1e4)
    S_unspot = Star(teff=T_unspot, radius=r_sun,mass=m_sun,R=1e4)
    
    s_spot = S_spot.spectrum(new.wavelength)/np.median(S_spot.spectrum(new.wavelength))
    s_unspot = S_unspot.spectrum(new.wavelength)/np.median(S_unspot.spectrum(new.wavelength))
    
    flux_ratio = s_spot/s_unspot
    
    top = (1. - f_tra) + f_tra * flux_ratio
    
    bottom = (1. - f_spot) + f_spot * flux_ratio
    
    delta_D_lam = (top / bottom) - 1
    
    depth = transit_depth * delta_D_lam
    
    return depth
