from ...imports import *
from ...spectra import *

__all__ = ["inject_spectrum"]


def inject_spectrum(
    self,
    temperature=5800 * u.K,
    logg=4.43,
    metallicity=0.0,
    radius=1 * u.Rsun,
    distance=10 * u.pc,
    phoenix=True,
):
    """
    Inject a stellar spectrum.

    Parameters
    ----------
    temperature : u.Quantity
        Temperature, in K (with no astropy units attached).
    logg : float
        Surface gravity log10[g/(cm/s**2)] (with no astropy units attached).
    metallicity : float
        Metallicity log10[metals/solar] (with no astropy units attached).
    radius : u.Quantity
        The radius of the star.
    distance : u.Quantity
        The distance to the star.
    phoenix : bool
        If `True`, use PHOENIX surface flux.
        If `False`, use Planck surface flux.

    Returns
    -------
    rainbow : Rainbow
        A new Rainbow object with the spectrum injected.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inject_spectrum", locals())

    # create a copy of the existing Rainbow
    new = self._create_copy()

    # warn if maybe we shouldn't inject anything
    if np.all(self.flux != 1):
        warnings.warn(
            f"""
        None of the pre-existing flux values were 1,
        which hints at the possibility that there
        might already be a spectrum in them. Please
        watch out for weird units or values!
        """
        )

    if phoenix:
        f = get_phoenix_photons
    else:
        f = get_planck_photons

    # get the spectrum from the surface
    _, surface_flux = f(
        temperature=u.Quantity(temperature).value,
        logg=logg,
        metallicity=metallicity,
        wavelength=self.wavelength,
    )

    # get the received flux at Earth
    received_flux = surface_flux * (radius / distance).decompose() ** 2

    # do math with spectrum
    for k in ["flux", "model", "uncertainty"]:
        try:
            new.fluxlike[k] = self.get(k) * self._broadcast_to_fluxlike(received_flux)
        except KeyError:
            pass

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new object
    return new
