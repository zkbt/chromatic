from ...imports import *

__all__ = [
    "inject_systematics",
    "_create_fake_wavelike_quantity",
    "_create_fake_timelike_quantity",
    "_create_fake_fluxlike_quantity",
]

"""
The idea for simple GP-generated systematics came from the fabulous
tutorial Neale Gibson led for the ers-transit team in summer 2021

The recording and a follow-along notebook are available at
https://ers-transit.github.io/pre-launch-hackathon.html
"""


def create_covariance_matrix(x, correlated=1, length=3, uncorrelated=0.5):
    """
    Generate a covariance matrix from a squared-exponential kernel
    with an additional diagonal component. The kernel is given by
    the following expression for x[m] and x[n]:

    distance = x[m] - x[n]
    correlated**2 * exp(-0.5*distance**2/length**2) + uncorrelated**2 * delta(m,n)

    It will be normalized to correlated**2 + uncorrelated**2 = 1

    Parameters
    ----------
    x : array or u.Quantity
        The arry of values for generating the matrix.
    correlated : float
        The amplitude of the correlated noise component.
    length : float or u.Quantity
        The length scale of the correlation (in same units as x)
    uncorrelated : float
        The amplitude of the uncorrelated noise component.
    """

    # what's the distance between each point and each other (an square matrix)
    distance_squared = np.subtract.outer(x, x) ** 2

    # calculate normalized coefficients
    N = correlated**2 + uncorrelated**2
    A, B = correlated / N, uncorrelated / N

    # generate the covariance matrix with squared exponential kernel
    covariance_matrix = A**2 * np.exp(-0.5 * distance_squared / length**2)

    # add extra variance along the diagonal
    i = np.arange(len(x))
    covariance_matrix[i, i] += B**2

    return covariance_matrix


def _create_fake_timelike_quantity(
    self, mean=np.zeros_like, correlated=1, dt=0.5 * u.hour, uncorrelated=0.25
):
    """
    Create a fake timelike quantity.

    Parameters
    ----------
    mean : function
        A function that takes in an array of times and returns
        a mean value for each time. The fake timeseries will
        wobble around that mean, according to the covariance
        matrix.
    correlated : float
        The amplitude of the correlated noise component.
    dt : Quantity
        The time scale of the correlation (in units of time).
    uncorrelated : float
        The amplitude of the uncorrelated noise component.
        (correlated and uncorrelated will be normalized to 1)
    Returns
    -------
    timelike : array
        A timelike array with the requested noise properties.
    """
    x = self.time
    covariance_matrix = create_covariance_matrix(
        x, correlated=correlated, length=dt, uncorrelated=uncorrelated
    )
    return np.random.multivariate_normal(mean(x), covariance_matrix)


def _create_fake_wavelike_quantity(
    self, mean=np.zeros_like, correlated=1, R=5, uncorrelated=0.25
):
    """
    Create a fake timelike quantity.

    Parameters
    ----------
    mean : function
        A function that takes in an array of times and returns
        a mean value for each time. The fake timeseries will
        wobble around that mean, according to the covariance
        matrix.
    correlated : float
        The amplitude of the correlated noise component.
    dt : Quantity
        The time scale of the correlation (in units of time).
    uncorrelated : float
        The amplitude of the uncorrelated noise component.
        (correlated and uncorrelated will be normalized to 1)

    Returns
    -------
    wavelike : array
        A wavelike array with the requested noise properties.
    """
    x = np.log(self.wavelength.value)
    covariance_matrix = create_covariance_matrix(
        x, correlated=correlated, length=1 / R, uncorrelated=uncorrelated
    )
    return np.random.multivariate_normal(mean(x), covariance_matrix)


def _create_fake_fluxlike_quantity(self, timelike_kw={}, wavelike_kw={}):
    """
    Create a fake fluxlike quantity.

    Parameters
    ----------
    timelike_kw : dict, optional
        Dictionary of keywords to pass to `_create_fake_timelike_quantity()`

    Returns
    -------
    fluxlike : array
        A fluxlike array (with )
    """

    # create a 1D time array
    tkw = dict(correlated=1, uncorrelated=0.5, dt=1 * u.hour)
    tkw.update(**timelike_kw)
    t = self._create_fake_timelike_quantity(**tkw)

    # create a 1D wavelength array
    wkw = dict(correlated=1, uncorrelated=0, R=5)
    wkw.update(**wavelike_kw)
    w = self._create_fake_wavelike_quantity(**wkw)

    # multiply them together and return
    f = w[:, np.newaxis] * t[np.newaxis, :]
    return f


def inject_systematics(
    self,
    amplitude=0.003,
    wavelike=[],
    timelike=["x", "y", "time"],
    fluxlike=["background"],
):

    """
    Inject some (very cartoony) instrumental systematics.

    Here's the basic procedure:

    1) Generate some fake variables that vary either just with
    wavelength, just with time, or with both time and wavelength.
    Store these variables for later use. For example, these might
    represent an average `x` and `y` centroid of the trace on the
    detector (one for each time), or the background flux associated
    with each wavelength (one for each time and for each wavelength).

    2) Generate a flux model as some function of those variables.
    In reality, we probably don't know the actual relationship
    between these inputs and the flux, but we can imagine one!

    3) Inject the model flux into the `flux` of this Rainbow,
    and store the combined model in `systematics-model` and
    each individual component in `systematic-model-{...}`.

    Parameters
    ----------
    amplitude : float, optional
        The (standard deviation-ish) amplitude of the systematics
        in units normalized to 1. For example, an amplitude of 0.003
        will produce systematic trends that tend to range (at 1 sigma)
        from 0.997 to 1.003.
    wavelike : list of strings, optional
        A list of wave-like cotrending quantities to serve as ingredients
        to a linear combination systematics model. Existing quantities
        will be pulled from the appropriate core dictionary; fake
        data will be created for quantities that don't already exist,
        from a cartoony Gaussian process model.
    timelike : list of strings, optional
        A list of time-like cotrending quantities to serve as ingredients
        to a linear combination systematics model. Existing quantities
        will be pulled from the appropriate core dictionary; fake
        data will be created for quantities that don't already exist,
        from a cartoony Gaussian process model.
    fluxlike : list of strings, optional
        A list of flux-like cotrending quantities to serve as ingredients
        to a linear combination systematics model. Existing quantities
        will be pulled from the appropriate core dictionary; fake
        data will be created for quantities that don't already exist,
        from a cartoony Gaussian process model.

    Returns
    -------
    rainbow : Rainbow
        A new Rainbow object with the systematics injected.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("inject_systematics", locals())

    # create a copy of the existing Rainbow
    new = self._create_copy()
    new.fluxlike["systematics_model"] = np.ones(self.shape)

    def standardize(q):
        """
        A quick helper to normalize all inputs to zero mean
        and unit standard deviation. It
        """
        offset = np.nanmean(q)
        sigma = np.nanstd(q)
        return u.Quantity((q - offset) / sigma).value, offset, sigma

    components = {}
    for k in wavelike:
        if k in self.wavelike:
            x, offset, sigma = standardize(self.wavelike[k])
        else:
            x = new._create_fake_wavelike_quantity()
            offset, sigma = 0, 1
            new.wavelike[k] = x
        c = np.random.normal(0, amplitude)
        df = c * x[:, np.newaxis] * np.ones(self.shape)
        new.fluxlike[f"systematics_model_from_{k}"] = df
        new.fluxlike["systematics_model"] += df
        components.update(
            **{
                f"linear_{k}": f"c_{k}*({k} - offset_{k})/sigma_{k}",
                f"c_{k}": c,
                f"offset_{k}": offset,
                f"sigma_{k}": sigma,
            }
        )

    for k in timelike:
        if k in self.timelike:
            x, offset, sigma = standardize(self.timelike[k])
        else:
            x = new._create_fake_timelike_quantity()
            offset, sigma = 0, 1
            new.timelike[k] = x
        c = np.random.normal(0, amplitude)
        df = c * x[np.newaxis, :] * np.ones(self.shape)
        new.fluxlike[f"systematics_model_from_{k}"] = df
        new.fluxlike["systematics_model"] += df
        components.update(
            **{
                f"linear_{k}": f"c_{k}*({k} - offset_{k})/sigma_{k}",
                f"c_{k}": c,
                f"offset_{k}": offset,
                f"sigma_{k}": sigma,
            }
        )

    for k in fluxlike:
        if k in self.fluxlike:
            x, offset, sigma = standardize(self.fluxlike[k])
        else:
            x = new._create_fake_fluxlike_quantity()
            offset, sigma = 0, 1
            new.fluxlike[k] = x
        c = np.random.normal(0, amplitude)
        df = c * x * np.ones(self.shape)
        new.fluxlike[f"systematics_model_from_{k}"] = df
        new.fluxlike["systematics_model"] += df
        components.update(
            **{
                f"linear_{k}": f"c_{k}*({k} - offset_{k})/sigma_{k}",
                f"c_{k}": c,
                f"offset_{k}": offset,
                f"sigma_{k}": sigma,
            }
        )

    new.metadata["systematics_components"] = components
    new.metadata["systematics_equation"] = "f = 1\n  + " + "\n  + ".join(
        [v for k, v in components.items() if k[:7] == "linear_"]
    )

    # modify both the model and flux arrays
    new.flux *= new.systematics_model
    new.model = new.fluxlike.get("model", 1) * new.systematics_model

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new object
    return new
