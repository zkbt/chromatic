from ...imports import *

__all__ = ["get_spectral_resolution"]


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
