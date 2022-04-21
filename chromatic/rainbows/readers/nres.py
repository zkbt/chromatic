"""
This module serves as a template for creating a new Rainbow
reader. If you want to add the ability to read chromatic
light curves from a new kind of file format, a good process
would be to do something like the following

    1. Copy this `template.py` file into a new file in the
    `readers/` directory, ideally with a name that's easy
    to recognize, such as `readers/abcdefgh.py` (assuming
    `abcdefgh` is the name of your format)

    2. Start by finding and replacing `abcdefgh` in this
    template with the name of your format.

    3. Edit the `from_abcdefgh` function so that it will
    load a chromatic light curve file in your format and,
    for some Rainbow object `rainbow`, populate at least:

        + rainbow.timelike['time']
        + rainbow.wavelike['wavelength']
        + rainbow.fluxlike['flux']

    You'll need to replace the cartoon functions on each
    line with the actual code needed to load your file.

    (This template assumes that only one file needs to be
    loaded. If you need to load multiple segments, or each
    time point is stored in its own file or something, then
    check out `stsci.py` for an example of loading and
    stitching together multiple input files. You'll probably
    want to change `filepath` to accept a glob-friendly string
    like `my-neato-formatted-files-*.npy` or some such.)

    4. Edit the `readers/__init__.py` file to import your
    `from_abcdefgh` function to be accessible when people
    are trying to create new Rainbows. Add an `elif` statement
    to the `guess_reader` function that will help guess which
    reader to use from some aspect(s) of the filename.

    (This `guess_reader` function also accepts a `format=`
    keyword that allows the user to explicitly specify that
    the abcdefgh reader should be used.)

    5. Submit a pull request to the github repository for
    this package, so that other folks can use your handy
    new reader too!
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_nres"]


def from_nres(rainbow, filepath):
    """
    Populate a Rainbow from a file in the NRES format.

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
            + `rainbow.timelike['time']`
            + `rainbow.wavelike['wavelength']`
            + `rainbow.fluxlike['flux']`
        and optionally, additional entries like
            + `rainbow.metadata['some-useful-parameter']`
            + `rainbow.fluxlike['uncertainty']`
            + `rainbow.fluxlike['ok']`

    filepath : str
        The path to the file to load.
        
    order : float
        The order to extract.
        Acceptable values span 52 to 119
    """

    # read in your file, however you like
    filenames = glob.glob(filepath)
    filenames = np.sort(filenames)
    
    hdu = fits.open(filenames[0]) # this will have to be updated
    
    order_index = order - 52 # the 0th order is order 52
    
    date = hdu['PRIMARY'].header['MJD-OBS']
    date_bjd = date + 2400000.5

    science_fiber_id = hdu['PRIMARY'].header['SCIFIBER']
    is_science = hdu['SPECTRUM'].data['FIBER'] == science_fiber_id
    data = hdu['SPECTRUM'].data[is_science]
    this_order = data[order_index]
    order = this_order['ORDER']
    mask = this_order['MASK'] == 0

    # these should be the good data for the particular order in question
    good_wavelengths = np.array([this_order['WAVELENGTH'][mask] / 10]) * u.nm
    good_normalized_fluxes = this_order['NORMFLUX'][mask]
    good_normalized_uncertainties = np.abs(this_order['NORMUNCERTAINTY'][mask])

    sort = np.argsort(good_wavelengths.value)

    wavelengths = good_wavelengths[sort]
    fluxes = good_normalized_fluxes[sort]
    errors = good_normalized_uncertainties[sort]
    

    # populate a 1D array of wavelengths (with astropy units of length)
    rainbow.wavelike["wavelength"] = wavelengths

    # populate a 1D array of times (with astropy units of time)
    rainbow.timelike["time"] = date_bjd

    # populate a 2D (row = wavelength, col = array of fluxes)
    rainbow.fluxlike["flux"] = np.array(wavelengths,fluxes)
    
    # populate a 2D (row = wavelength, col = array of errors)
    rainbow.fluxlike["uncertainty"] = np.array(wavelengths, errors)
    
    rainbow.fluxlike["mask"] = mask

    # add some warnings if there's any funny business
    if len(filenames) == 0:
        warnings.warn(
            f"""
        There are no files of that name in this folder!
        """
        )