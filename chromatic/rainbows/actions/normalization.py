from ...imports import *


def normalize(self, wavelength=True, time=False):
    """
    Normalize by dividing through by the median spectrum and/or lightcurve.

    Parameters
    ----------
    wavelength : bool
        Should we divide by the median spectrum?

    time : bool
        Should we divide by the median light curve?

    Returns
    -------
    normalized : MultiRainbow
        The normalized MultiRainbow.
    """

    # TODO, think about more careful treatment of uncertainties + good/bad data
    new = self._create_copy()

    # (ignore nan warnings)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if time:
            new = new / np.nanmedian(self.flux, axis=self.waveaxis)

        if wavelength:
            new = new / np.nanmedian(self.flux, axis=self.timeaxis)

    return new
