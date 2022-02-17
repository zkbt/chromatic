"""
Basic model injection capability. *Very* in-progress.
"""
from ..imports import *
import starry


def make_fake_radii(mean=0.1, wobble=0.005):
    def kludge_radii(wavelength):
        phase = wavelength.to("micron").value
        return mean + wobble * np.sin(2 * np.pi * phase)

    return kludge_radii


def make_fake_limb_darkening(max_wavelength=5 * u.micron):
    def kludge_limb_darkening(wavelength):
        u1 = (max_wavelength - wavelength) / max_wavelength
        u2 = 0  # 0.5*(max_wavelength - wavelength)/max_wavelength
        return np.minimum(np.maximum(np.vstack([u1, u2]).value, 0), 1)

    return kludge_limb_darkening


def make_fake_flux_ratio(T_planet=1500 * u.K, T_star=5000 * u.K):
    import rainbowconnection as rc

    def kludge_flux_ratio(wavelength):
        B_planet = rc.Thermal(T_planet)
        B_star = rc.Thermal(T_star)
        return (
            B_planet.surface_flux(wavelength) / B_star.surface_flux(wavelength)
        ).value

    return kludge_flux_ratio


class StarryModelInjector:
    """
    This should be used only for injecting transit/phase-curve/eclipse
    models. It's not necessarily set up for easy optimization
    of the model parameters through comparison to data.
    """

    def __init__(
        self,
        planet_radius=make_fake_radii(),
        planet_flux_ratio=make_fake_flux_ratio(),
        limb_darkening=make_fake_limb_darkening(),
        stellar_mass=1.0,
        stellar_radius=1.0,
        inc=90.0,
        porb=1.0,
        t0=0.0,
        ecc=0.0,
        w=90.0,
        phase_curve_amplitude=0.5,
        phase_curve_offset=0,
    ):

        """
        Given some planet and stellar parameters, set up a
        starry model that will be stored inside this object.
        Once the model is defined, the model can be
        evaluated for a particular set of times/wavelengths.

        FIXME:
        + add Parameters docs
        + add better handling for wavelength-dependent quantities
            + allow function or array or scalar?
        + create an inference version of this object
            + set up pmcy3 model
            + test performance of one System per wavelength
            + build interface for fitted vs fixed
                + shared, wavelike, (timelike?)

        """
        starry.config.lazy = False

        # create a simple star
        pristine_star_map = starry.Map(
            ydeg=0,  # no asymmetries
            udeg=2,  # quadratic limb-darkening
            amp=1.0,  # set brightness scale to 1
        )
        self.star = starry.Primary(
            pristine_star_map, m=stellar_mass, r=stellar_radius, prot=1.0
        )

        # create a simple planet
        planet_map = starry.Map(
            ydeg=1,  # dipole map for the planet
            udeg=0,  # no limb darkening for planet
            amp=0.001,  # Fp/F*
            inc=90.0,  # edge-on inclination (is this relative to orbit, or us?)
            obl=0.0,
        )  # no obliquity
        self.planet = starry.Secondary(
            planet_map,
            m=0.0,  # ignore mass of planet
            r=0.1,  # rp/rs
            inc=inc,  # edge-on inclination
            prot=1.0,  # rotation of planet (is this days, or relative to period?)
            porb=porb,  # orbital period
            t0=t0,  # time of a reference transit
            ecc=ecc,  # eccentricity
            w=w,  # argument of periastron, in degrees
            theta0=180 + phase_curve_offset,
        )
        self.planet.map[1, 0] = phase_curve_amplitude

        # put the star and planet together into a system
        self.system = starry.System(self.star, self.planet)

        # set the wavelength-dependent functions
        self.planet_radius = planet_radius
        self.planet_flux_ratio = planet_flux_ratio
        self.limb_darkening = limb_darkening

    def __call__(self, rainbow):
        """
        When calling this model with a rainbow object,
        return another matched rainbow populated with model fluxes.
        """

        # create a copy of the rainbow
        new = rainbow._create_copy()

        # loop through wavelengts to populate fluxes
        t = new.time.to("day").value
        for i, w in enumerate(new.wavelength):
            self.system.secondaries[0].r = self.planet_radius(w)
            self.system.secondaries[0].map.amp = (
                self.planet_flux_ratio(w) * self.system.secondaries[0].r ** 2
            )

            u1, u2 = self.limb_darkening(w)
            self.system.primary.map[1] = u1
            self.system.primary.map[2] = u2
            new.flux[i, :] = self.system.flux(t)

        return new
