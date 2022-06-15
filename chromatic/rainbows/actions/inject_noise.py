from ...imports import *

__all__ = ["inject_noise"]


def inject_noise(self, signal_to_noise=100):
    """
    Inject uncorrelated Gaussian random noise into the flux.

    Parameters
    ----------
    signal_to_noise : float
        The signal-to-noise per wavelength per time.
        For example, S/N=100 would mean that the
        uncertainty on the flux for each each
        wavelength-time data point will be 1%.

    Returns
    -------
    rainbow : Rainbow
        A new Rainbow object with the systematics injected.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inject_noise", locals())

    # create a copy of the existing Rainbow
    new = self._create_copy()

    # get the underlying model (or create one if needed)
    if "model" in self.fluxlike:
        model = self.fluxlike["model"]
    else:
        # kludge, do we really want to allow this?
        model = self.flux * 1
        self.fluxlike["model"] = model

    # calculate the uncertainty with a fixed S/N
    uncertainty = model / signal_to_noise
    new.fluxlike["uncertainty"] = uncertainty

    # inject a realization of the noise
    new.fluxlike["flux"] = np.random.normal(model, uncertainty)

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new object
    return new
