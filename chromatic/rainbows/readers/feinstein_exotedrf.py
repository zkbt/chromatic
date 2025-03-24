# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_feinstein_exotedrf"]


def from_feinstein_exotedrf(rainbow, filepath):
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
    d = np.load(filepath, allow_pickle=True)[()]

    rainbow.wavelike["wavelength"] = d['wavelength'] * u.micron * 1

    # populate a 1D array of times (with astropy units of time)
    rainbow.timelike["time"] =  d['time'] * u.day * 1

    # populate a 2D (row = wavelength, col = array of fluxes
    flux = np.zeros((len(d['wavelength']), len(d['time'])))
    uncertainty = np.zeros_like(flux)

    rainbow.fluxlike["flux"] = d['flux'].T * 1
    rainbow.fluxlike["uncertainty"] = d['flux_err'].T * 1
