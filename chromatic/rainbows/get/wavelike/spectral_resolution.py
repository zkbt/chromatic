from ....imports import *

__all__ = ["get_spectral_resolution"]


def get_spectral_resolution(self, pixels_per_resolution_element=1):
    """
    Estimate the R=w/dw spectral resolution.

    Higher spectral resolutions correspond to more wavelength
    points within a particular interval. By default, it's
    estimated for the interval between adjacent wavelength
    bins. In unbinned data coming directly from a telescope,
    there's a good chance that adjacent pixels both sample
    the same resolution element as blurred by the telescope
    optics, so the `pixels_per_resolution_element` keyword
    should likely be larger than 1.

    Parameters
    ----------
    pixels_per_resolution_element : float, optional
        How many pixels do we consider as a resolution element?

    Returns
    -------
    R : array
        The spectral resolution at each wavelength.
    """

    # calculate spectral resolution, for this pixels/element
    w = self.wavelength
    dw = np.gradient(self.wavelength)
    R = np.abs(w / dw / pixels_per_resolution_element)

    return R
