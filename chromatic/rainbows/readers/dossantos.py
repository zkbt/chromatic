"""
This module serves as a template for creating a new Rainbow
reader. If you want to add the ability to read chromatic
light curves from a new kind of file format, a good process
would be to do something like the following

    1. Copy this `template.py` file into a new file in the
    `readers/` directory, ideally with a name that's easy
    to recognize, such as `readers/dossantos.py` (assuming
    `dossantos` is the name of your format)

    2. Start by finding and replacing `dossantos` in this
    template with the name of your format.

    3. Edit the `from_dossantos` function so that it will
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
    `from_dossantos` function to be accessible when people
    are trying to create new Rainbows. Add an `elif` statement
    to the `guess_reader` function that will help guess which
    reader to use from some aspect(s) of the filename.

    (This `guess_reader` function also accepts a `format=`
    keyword that allows the user to explicitly specify that
    the dossantos reader should be used.)

    5. Submit a pull request to the github repository for
    this package, so that other folks can use your handy
    new reader too!
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_dossantos"]


def from_dossantos(rainbow, filepath):
    """
    Populate a Rainbow from a file in the dossantos format.

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
    """

    import pickle

    with open(filepath, "rb") as pkl:
        spectra = pickle.load(pkl)

    # populate a 1D array of wavelengths (with astropy units of length)
    wavelength = spectra[0]["wavelength"] * u.micron
    rainbow.wavelike["wavelength"] = wavelength * 1

    # populate a 1D array of times (with astropy units of time)
    times = np.arange(len(spectra)) * u.minute
    cheerfully_suggest("The times are totally made up!")
    rainbow.timelike["time"] = times * 1

    # populate a 2D (row = wavelength, col = array of fluxes
    flux = np.zeros((len(wavelength), len(times)))
    uncertainty = np.zeros_like(flux)

    for k in spectra:
        flux[:, k] = spectra[k]["flux"]
        uncertainty[:, k] = spectra[k]["flux_error"]

    rainbow.fluxlike["flux"] = flux * 1
    rainbow.fluxlike["uncertainty"] = uncertainty * 1
