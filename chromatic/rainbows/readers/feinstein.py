# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_feinstein_numpy", "from_feinstein_h5"]


def from_feinstein_numpy(rainbow, filepath):
    """
    Populate a Rainbow from a file in the feinstein numpy format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load.
    """

    try:
        from h5py import File
    except ImportError:
        warnings.warn(
            f"""
        Please try to install `h5py` into your current environment with
            `pip install --upgrade h5py`
        to have access to the `h5py` reader needed for HDF5 files.
        """
        )
    time, wavelength, spectra, err = np.load(filepath, allow_pickle=True)

    rainbow.wavelike["wavelength"] = wavelength * u.micron * 1

    # populate a 1D array of times (with astropy units of time)
    times = time * u.day
    rainbow.timelike["time"] = times * 1

    # populate a 2D (row = wavelength, col = array of fluxes
    flux = np.zeros((len(wavelength), len(times)))
    uncertainty = np.zeros_like(flux)

    rainbow.fluxlike["flux"] = spectra.T * 1
    rainbow.fluxlike["uncertainty"] = err.T * 1


def from_feinstein_h5(self, filepath, order=1, version="opt"):
    """
    Populate a Rainbow from a file in the feinstein H5 format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load.
    """

    f = File(filepath)

    astropy_times = Time(np.array(f["time"]), format="mjd", scale="tdb")
    self.set_times_from_astropy(astropy_times, is_barycentric=True)  # ???
    self.wavelike["wavelength"] = (
        np.array(f[f"wavelength_order_{order}"]) * u.micron * 1
    )
    for k in ["box_flux", "box_var", "opt_flux", "opt_var"]:
        self.fluxlike[k] = np.transpose(np.array(f[f"{k}_order_{order}"]) * 1)

    self.fluxlike["flux"] = self.fluxlike[f"{version}_flux"] * 1
    self.fluxlike["uncertainty"] = self.fluxlike[f"{version}_var"] * 1
