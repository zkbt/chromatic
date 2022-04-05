from ...imports import *


def normalize(self, wavelength=True, time=False, percentile=50):
    """
    Normalize by dividing through by the median spectrum and/or lightcurve.

    Parameters
    ----------
    wavelength : bool
        Should we divide by the median spectrum?

    time : bool
        Should we divide by the median light curve?

    percentile : float
        A number between 0 and 100, specifying the percentile
        of the data along an axis to use as the reference.
        The default of `percentile=50` corresponds to the median.
        If you want to normalize to out-of-transit, maybe you
        want a higher percentile. If you want to normalize to
        the baseline below a flare, maybe you want a lower
        percentage.

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
            new = new / np.nanpercentile(self.flux, axis=self.waveaxis)

        if wavelength:
            new = new / np.nanpercentile(self.flux, axis=self.timeaxis)

    return new
