# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_radica"]


def from_radica(self, filepath, order=1):
    """
    Populate a Rainbow from a file in the radica format.

    Parameters
    ----------

    self : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load.
    """

    hdu = fits.open(filepath)

    n_orders = 2
    cheerfully_suggest(
        f"""
    Loading NIRISS spectroscopic `order={order}``. Two orders are available,
    and you can set which (1,2) you want to read with the `order=` option.
    """
    )
    self.fluxlike["wavelength_2d"] = (
        np.transpose(hdu[f"Wave 2D Order {order}"].data * u.micron) * 1
    )
    self.wavelike["wavelength"] = (
        np.nanmedian(self.fluxlike["wavelength_2d"], axis=1) * 1
    )
    self.fluxlike["flux"] = (
        np.transpose(hdu[f"Flux Order {order}"].data * u.electron) * 1
    )
    self.fluxlike["uncertainty"] = (
        np.transpose(hdu[f"Flux Error Order {order}"].data * u.electron) * 1
    )
    astropy_times = Time(hdu["Time"].data, format="jd", scale="tdb")
    self.set_times_from_astropy(astropy_times, is_barycentric=True)
