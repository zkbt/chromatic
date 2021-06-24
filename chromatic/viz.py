#Code to create plots of rainbows

import matplotlib.pyplot as plt
import numpy as np


def wavelength_plot(flux, time, wavelength, step_size = None, **kwargs):
    """
    Plot a time series for each wavelength bin

    Parameters
    ----------
    flux : TYPE
        DESCRIPTION.
    time : TYPE
        DESCRIPTION.
    wavelength : TYPE
        DESCRIPTION.
    step_size : float, default=None
        Spacing between each time series. None defaults to half the standard dev.

    Returns
    -------
    None.

    """
    
    
    mean_flux = np.nanmean(flux, axis=0)
    median_flux = np.nanmedian(flux, axis=0)
    
    min_flux = np.min(median_flux)
    max_flux = np.max(median_flux)
    
    if step_size is None:
        step_size = 0.5 * np.std(flux)
    
    
    nsteps = len(wavelength)
    
    for i, w in enumerate(wavelength):
        
        lc = flux[i, :]
        plot_flux = lc
        
        color = (0, 0.2 - 0.2 * (i / nsteps), i / nsteps)
        plt.figure(figsize = (10,8))
        
        plt.plot(time, plot_flux, '.', color = color)
        