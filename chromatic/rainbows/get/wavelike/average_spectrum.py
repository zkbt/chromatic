from ....imports import *

__all__ = ["get_average_spectrum"]


def get_average_spectrum(self):
    """
    Return a average_spectrum of the star, averaged over all times.

    This uses `bin`, which is a horribly slow way of doing what is
    fundamentally a very simply array calculation, because we
    don't need to deal with partial pixels.

    Returns
    -------
    average_spectrum : array
        Wavelike array of average spectrum.
    """
    return self.get_average_spectrum_as_rainbow().flux[:, 0]
