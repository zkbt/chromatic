"""
Read multiple exposures of a single order of LCO/NRES spectra.
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_nres"]


def from_nres(rainbow, filepath, order=52):
    """
    Populate a Rainbow from a file in the NRES format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.

    filepath : str
        The path to the file to load.

    order : float
        The order to extract.
        Acceptable values span 52 to 119
    """

    # read in your file, however you like
    filenames = glob.glob(filepath)
    filenames = np.sort(filenames)

    order_index = order - 52  # the 0th order is order 52

    for i, f in tqdm(enumerate(filenames), leave=False):

        # open this
        hdu = fits.open(f)

        # get the time associated with this observation
        date = hdu["PRIMARY"].header["MJD-OBS"]
        date_bjd = date + 2400000.5

        # get the science spectrum
        science_fiber_id = hdu["PRIMARY"].header["SCIFIBER"]
        is_science = hdu["SPECTRUM"].data["FIBER"] == science_fiber_id
        data = hdu["SPECTRUM"].data[is_science]
        this_order = data[order_index]
        order = this_order["ORDER"]
        ok = this_order["MASK"] == 0

        # these should be the good data for the particular order in question
        wavelengths = np.array(this_order["WAVELENGTH"]) * u.angstrom  # [mask]
        normalized_fluxes = this_order["NORMFLUX"]  # [mask]
        normalized_uncertainties = np.abs(this_order["NORMUNCERTAINTY"])  # [mask])

        sort = np.argsort(wavelengths.value)
        wavelengths = wavelengths[sort]
        fluxes = normalized_fluxes[sort]
        errors = normalized_uncertainties[sort]

        if i == 0:
            ntimes = len(filenames)
            nwaves = len(wavelengths)
            for k in ["wavelength_2d", "flux", "uncertainty", "ok"]:
                rainbow.fluxlike[k] = np.zeros((nwaves, ntimes))
            rainbow.fluxlike["wavelength_2d"] *= u.micron
            rainbow.fluxlike["ok"] = rainbow.fluxlike["ok"].astype(np.bool)
            rainbow.timelike["time"] = np.zeros(ntimes) * u.day

        # populate a 1D array of times (with astropy units of time)
        rainbow.timelike["time"][i] = date_bjd * u.day * 1

        # populate a 2D (row = wavelength, col = time, value = wavelength)
        rainbow.fluxlike["wavelength_2d"][:, i] = wavelengths * 1

        # populate a 2D (row = wavelength, col = time) array of fluxes
        rainbow.fluxlike["flux"][:, i] = fluxes * 1

        # populate a 2D (row = wavelength, col = time) array of uncertainties
        rainbow.fluxlike["uncertainty"][:, i] = errors * 1

        # populate a 2D (row = wavelength, col = time) array of ok
        rainbow.fluxlike["ok"][:, i] = ok * 1

    # populate a 1D array of wavelengths (with astropy units of length)
    rainbow.wavelike["wavelength"] = np.nanmedian(
        rainbow.fluxlike["wavelength_2d"], axis=rainbow.timeaxis
    )

    # add some warnings if there's any funny business
    if len(filenames) == 0:
        cheerfully_suggest(
            f"""
        There are no files of that name in this folder!
        """
        )
