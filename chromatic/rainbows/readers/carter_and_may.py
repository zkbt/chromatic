# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_carter_and_may"]


def from_carter_and_may(self, filepath):
    """
    Populate a Rainbow from a file in the ERS benchmark
    transmission spectrum of WASP-39b format.

    https://zenodo.org/records/10161743

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated. This function is meant
        to be called during the initialization of a Rainbow
        object. You can/should assume that the `rainbow` object
        being passed will already contain the four core
        dictionaries `timelike`, `wavelike`, `fluxlike`, and
        `metadata`. This function should, at a minimum, add
        the following items
            + `self.timelike['time']`
            + `self.wavelike['wavelength']`
            + `self.fluxlike['flux']`
        and optionally, additional entries like
            + `self.metadata['some-useful-parameter']`
            + `self.fluxlike['uncertainty']`
            + `self.fluxlike['ok']`

    filepath : str
        The path to the file to load.
    """

    # read in your file, however you like
    f = h5.File(filepath)

    # populate a 1D array of wavelengths (with astropy units of length)
    self.wavelike["wavelength"] = f["wave_1d"] * u.micron

    # populate a 1D array of times (with astropy units of time)
    self.timelike["time"] = (f["time"]) * u.day + 2400000.5 * u.day

    # populate a 2D (row = wavelength, col = time) array of fluxes
    self.fluxlike["flux"] = np.transpose(f["optspec"])
    self.fluxlike["uncertainty"] = np.transpose(f["opterr"])
    self.fluxlike["ok"] = np.transpose(f["optmask"]) == 0

    for k in f.keys():
        if k not in ["wave_1d", "time", "optspec", "opterr", "optmask"]:
            if np.shape(f[k]) == np.shape(self.time):
                self.timelike[k] = f[k]
            elif np.shape(f[k]) == np.shape(self.wavelength):
                self.wavelike[k] = f[k]
            elif np.shape(f[k]) == np.shape(self.flux.T):
                self.fluxlike[k] = np.transpose(f[k])

    ok = self.wavelength != 0
    for k in self.wavelike:
        self.wavelike[k] = self.wavelike[k][ok]
    for k in self.fluxlike:
        self.fluxlike[k] = self.fluxlike[k][ok, :]
