from ...imports import *
import batman

__all__ = ["inject_transit"]


def inject_transit(self, planet_params={}, planet_radius=0.1):

    """
    Simulate a wavelength-dependent planetary transit using
    batman.

    Parameters
    ----------

    planet_params : Dictionary
        Values for planetary parameters for use in batman modelling.
        Any values not supplied will be set to defaults:
            "t0" = time of inferior conjunction (days) (default 0)
            "per" = orbital period (days) (detault 1)
            "a" = semi-major axis (units of stellar radii) (default 15)
            "inc" = inclination (degrees) (default 90)
            "ecc" = eccentricity (default 0)
            "w" = longitude of periastron (degrees)(default 0)
            "limb_dark" = limb-darkening model (default "nonlinear"), possible
                values described in more detail in batman documentation
            "u" = limb-darkening coefficients (default [0.5, 0.1, 0.1, -0.1])
                Can take 3 forms:
                -A single value (if limb-darkening law requires only one value)
                -A 1D list/array of coefficients corresponding to the limb-darkening
                law
                -A 2D array of the form (n_wavelengths, n_coefficients) where
                each row is the set of limb-darkening coefficients corresponding
                to a single wavelength

                Note that this currently does not calculate the appropriate
                coefficient vs wavelength variations itself- there exist codes
                (such as hpparvi/PyLDTk and nespinoza/limb-darkening) which
                can be used for this.

        example value: planet_params = {"a":12, "inc":87}

    planet_radius = Two options:
            1D array with same dimensions as wavelength array,
                each value corresponds to planet radius/stellar radius at that
                wavelength.
            float representing Rp/Rstar if the radius is not wavelength-dependent.
        example value: planet_radius = 0.01,

    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inject_transit", locals())

    # create a copy of the existing Rainbow
    new = self._create_copy()

    # First, make sure planet_radius has the right dimension.
    if type(planet_radius) != float and len(planet_radius) != self.nwave:
        print(
            "Invalid planet radius array: must be float or have shape "
            + str(np.shape(self.wavelike["wavelength"]))
        )

    # Defaults for planet simulation.
    defaults = {
        "t0": 0,
        "per": 3,
        "a": 10,
        "inc": 90,
        "ecc": 0,
        "w": 0,
        "limb_dark": "nonlinear",
        "u": [0.5, 0.1, 0.1, -0.1],
    }

    # Read in planet parameters.
    for i in range(len(planet_params.keys())):
        key = list(planet_params.keys())[i]
        if key in list(defaults.keys()):
            defaults[key] = planet_params[key]
        else:
            print("Warning: " + str(key) + " not a valid parameter")

    # Initialize batman model.
    params = batman.TransitParams()
    params.t0 = defaults["t0"]
    params.per = defaults["per"]
    params.a = defaults["a"]
    params.inc = defaults["inc"]
    params.ecc = defaults["ecc"]
    params.w = defaults["w"]
    params.limb_dark = defaults["limb_dark"]

    # Deal with limb-darkening.
    if len(np.shape(defaults["u"])) < 2:  # Coefficients constant with wavelength
        u_arr = np.tile(defaults["u"], (self.nwave, 1))

    elif (
        len(np.shape(defaults["u"])) == 2
    ):  # 2D array of coefficients, along wavelength axis
        if np.shape(defaults["u"])[0] != self.nwave:
            print("Shape of limb-darkening array does not match wavelengths.")
            return
        u_arr = defaults["u"]

    else:
        print("Invalid limb-darkening coefficient array.")

    # Read in planetary radius.
    if type(planet_radius) == float:
        rprs = np.zeros(self.nwave) + planet_radius
    else:
        rprs = planet_radius

    planet_flux = np.zeros((self.nwave, self.ntime))

    for i in range(self.nwave):
        params.rp = rprs[i]
        params.u = u_arr[i]
        # print(params.u)
        try:
            m
        except NameError:
            m = batman.TransitModel(params, new.time.to_value("day"))
        planet_flux[i] = m.light_curve(params)

    new.planet_model = planet_flux
    new.flux *= new.planet_model
    new.model = new.fluxlike.get("model", 1) * new.planet_model

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new object
    return new
