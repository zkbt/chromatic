from ...imports import *

__all__ = ["inject_noise"]


def inject_noise(self, signal_to_noise=100, number_of_photons=None):
    """
    Inject uncorrelated random noise from a Gaussian
    or Poisson distribution into the flux.

    Parameters
    ----------

    signal_to_noise : float
        The signal-to-noise per wavelength per time.
        For example, S/N=100 would mean that the
        uncertainty on the flux for each each
        wavelength-time data point will be 1%.

    number_of_photons : float
        The number of photons expected to be recieved
        from some light source.  For example, the number
        of photons expected from a Sun-like star as seen
        with a human eye from a distance of one parsec and
        an exposure time of 0.1s is 25229 photons.

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
    if "model" in new.fluxlike:
        model = new.fluxlike["model"]
    else:
        # kludge, do we really want to allow this?
        model = self.flux * 1
        new.fluxlike["model"] = model

    # setting up an if/else statement so that the user
    # can choose if they want to use their own
    # number_of_photons or the automatic signal_to_noise
    # noise generation
    if number_of_photons != None:
        mu = model * number_of_photons

        # inject a realization of noise using number_of_photons
        # (yields poisson distribution)
        new.fluxlike["flux"] = np.random.poisson(mu) * u.photon #mu is the center

        # store S/N as metadata
        new.metadata["number_of_photons"] = number_of_photons

        # append the history entry to the new Rainbow
        new._record_history_entry(h)

    else:
        # calculate the uncertainty with a fixed S/N
        uncertainty = model / signal_to_noise
        new.fluxlike["uncertainty"] = uncertainty

        # inject a realization of the noise
        new.fluxlike["flux"] = np.random.normal(model, uncertainty)

        # store S/N as metadata
        new.metadata["signal_to_noise"] = signal_to_noise

        # append the history entry to the new Rainbow
        new._record_history_entry(h)

    # return the new object
    return new
