from ...imports import *

__all__ = ["inject_noise"]

# Molly! add a `number_of_photons=None` keyword argument to the function call.
# By saying, `=None`, it will default to None.
def inject_noise(self, signal_to_noise=100, number_of_photons=None):
    """
    Inject uncorrelated Gaussian random noise into the flux.

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
    print('Run ralphie run')
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

    # Molly! The next three chunks of code are the thing that would
    # need to be different for `number_of_photons`. I think you should
    # duplicate these to start with, and then set up an if/else statement
    # that does your `number_of_photons` calculation if its value
    # is not None, or does the `signal_to_noise` calculation otherwise.

    if number_of_photons != None:
        #mu = np.mean(number_of_photons)
        #std = np.sqrt(mu)
        mu = model * number_of_photons

        #uncertainty =
        #new.fluxlike["uncertainty"] = uncertainty
        #new.fluxlike["flux"] = np.random.poisson(mu, 100) #center, number of points)
        new.fluxlike["flux"] = np.random.poisson(mu)

        plt.hist(new.fluxlike["flux"].flatten(), bins = 30, density = True)

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
