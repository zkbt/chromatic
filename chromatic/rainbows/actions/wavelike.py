from ...imports import *


def get_spectral_resolution(self, pixels_per_resolution_element=1):
    """
    Parameters
    ----------
    pixels_per_resolution_element : float
        How many pixels do we consider as a resolution element?

    Returns
    -------
    R : np.array (wavelike)
        The spectral resolution, at each wavelength.
    """

    # calculate spectral resolution, for this pixels/element
    w = self.wavelength
    dw = np.gradient(self.wavelength)
    R = np.abs(w / dw / pixels_per_resolution_element)

    return R


def plot_spectral_resolution(
    self,
    pixels_per_resolution_element=1,
    ax=None,
    w_unit="micron",
    cmap=None,
    vmin=None,
    vmax=None,
    plotkw={},
    **kw,
):
    """
    Plot the spectral resolution as a function of wavelength.

    Parameters
    ----------
    pixels_per_resolution_element : float
        How many pixels do we consider as a resolution element?
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    w_unit : str, astropy.unit.Unit
        The unit for plotting wavelengths.
    cmap : str, matplotlib.colors.Colormap
        The color map to use for expressing wavelength.
    vmin : astropy.units.Quantity
        The minimum value to use for the wavelength colormap.
    vmax : astropy.units.Quantity
        The maximum value to use for the wavelength colormap.
    plotkw : dict
        A dictionary of keywords passed to `plt.plot`
    """
    R = self.get_spectral_resolution(pixels_per_resolution_element)
    w_unit = u.Unit(w_unit)

    with quantity_support():

        # make sure that the wavelength-based colormap is defined
        self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)

        # make sure ax is set up
        if ax is None:
            ax = plt.subplot()
        plt.sca(ax)

        # define default visuals, but let user override with scatterkw
        # this_scatterkw = dict(c=self.get_wavelength_color(self.wavelength))
        # this_scatterkw.update(**scatterkw)

        plt.plot(self.wavelength.to(w_unit), R, **plotkw)
        plt.xlabel(f'Wavelength ({w_unit.to_string("latex_inline")})')
        plt.ylabel(f"$R=\lambda/d\lambda$ ({pixels_per_resolution_element} pixel)")
        return ax

        ax.set_xlabel(f"Time ({self.time.unit.to_string('latex_inline')})")


def get_typical_uncertainty(self, function=np.nanmedian):
    """
    Get the typical per-wavelength uncertainty.

    Parameters
    ----------
    function : function
        What function should be used to choose the "typical"
        value for each wavelength? Good options are probably
        things like `np.nanmedian`, `np.median`, `np.nanmean`
        `np.mean`

    Returns
    -------
    uncertainty_per_wavelength : np.array (wavelike)
        The uncertainty associated with each wavelength.
    """
    uncertainty_per_wavelength = function(self.uncertainty, axis=self.timeaxis)
    return uncertainty_per_wavelength


def get_spectrum(self):
    """
    Return a spectrum of the star, averaged over all times.

    Parameters
    ----------


    Returns
    -------
    s : np.array (wavelike)
    """

    # FIXME, think about error weighting/masking
    return np.nanmedian(self.flux, axis=self.timeaxis)
