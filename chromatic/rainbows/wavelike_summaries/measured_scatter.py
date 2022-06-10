from ...imports import *

__all__ = ["get_measured_scatter"]


def get_measured_scatter(self, method='standard-deviation'):
    """
    Get measured scatter from a chromatic rainbow.
    
    Parameters
    ----------
    method : string
        What method to use to obtain measured scatter. Current options are 'MAD', 'standard-deviation'.
    
    Returns
    -------
    scatter : np.array (wavelike)
        The scatter wavelike array.
    """
    
    if method == 'standard-deviation':
        scatters = np.zeros(len(rainbow.wavelength))
        for i, row in enumerate(rainbow.flux):
            scatters[i-1] = np.std(row)
        return scatters
    
    if method == 'MAD':
        scatters = median_abs_deviation(rainbow.flux.copy(), axis=1, scale = 'normal') # Setting the scale 'normal' is like * 1.48
        return scatters
    
    else:
        print("Not an available method, choose from 'MAD', 'standard-deviation'.")