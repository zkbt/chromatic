from ...imports import *

__all__ = ["inject_noise"]


def inject_noise(self, signal_to_noise=100, number_of_photons=None):
    """
    Inject uncorrelated random noise into the `.flux` array.

    This injects independent noise to each data point,
    drawn from either a Gaussian or Poisson distribution.
    If the inputs can be scalar, or they can be arrays
    that we will try to broadcast into the shape of the
    `.flux` array.

    Parameters
    ----------

    signal_to_noise : float, array, optional
        The signal-to-noise per wavelength per time.
        For example, S/N=100 would mean that the
        uncertainty on the flux for each each
        wavelength-time data point will be 1%.
        If it is a scalar, then even point is the same.
        If it is an array with a fluxlike, wavelike,
        or timelike shape it will be broadcast
        appropriately.
    number_of_photons : float, array, optional
        The number of photons expected to be recieved
        from the light source per wavelength and time.
        If it is a scalar, then even point is the same.
        If it is an array with a fluxlike, wavelike,
        or timelike shape it will be broadcast
        appropriately.
        If `number_of_photons` is set, then `signal_to_noise`
        will be ignored.

    Returns
    -------
    rainbow : Rainbow
        A new `Rainbow` object with the noise injected.
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
    if number_of_photons is not None:
        if u.Quantity(model).unit != u.Unit(""):
            raise ValueError(
                f"""
            We haven't yet implemented `number_of_photons` noise
            for models that have units associated with them. Sorry!
            """
            )

        mu = model * self._broadcast_to_fluxlike(number_of_photons)

        # convert the model to photons and store it
        new.fluxlike["model"] = mu * u.photon

        # inject a realization of noise using number_of_photons
        # (yields poisson distribution)
        new.fluxlike["flux"] = np.random.poisson(mu) * u.photon  # mu is the center

        # store number of photons as metadata
        new.metadata["number_of_photons"] = number_of_photons

        # calculate the uncertainty
        uncertainty = np.sqrt(mu)
        new.fluxlike["uncertainty"] = uncertainty * u.photon

        # append the history entry to the new Rainbow
        new._record_history_entry(h)

    else:
        # calculate the uncertainty with a fixed S/N
        uncertainty = model / self._broadcast_to_fluxlike(signal_to_noise)
        new.fluxlike["uncertainty"] = uncertainty

        # inject a realization of the noise
        if isinstance(model, u.Quantity):
            unit = model.unit
            loc = model.to_value(unit)
            scale = uncertainty.to_value(unit)
        else:
            unit = 1
            loc = model
            scale = uncertainty
        new.fluxlike["flux"] = np.random.normal(model, uncertainty) * unit

        # store S/N as metadata
        new.metadata["signal_to_noise"] = signal_to_noise

        # append the history entry to the new Rainbow
        new._record_history_entry(h)

    # return the new object
    return new
