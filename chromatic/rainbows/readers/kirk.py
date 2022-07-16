"""
This module serves as a template for creating a new Rainbow
reader. If you want to add the ability to read chromatic
light curves from a new kind of file format, a good process
would be to do something like the following

    1. Copy this `template.py` file into a new file in the
    `readers/` directory, ideally with a name that's easy
    to recognize, such as `readers/kirk.py` (assuming
    `kirk` is the name of your format)

    2. Start by finding and replacing `kirk` in this
    template with the name of your format.

    3. Edit the `from_kirk` function so that it will
    load a chromatic light curve file in your format and,
    for some Rainbow object `rainbow`, populate at least:

        + self.timelike['time']
        + self.wavelike['wavelength']
        + self.fluxlike['flux']

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
    `from_kirk` function to be accessible when people
    are trying to create new Rainbows. Add an `elif` statement
    to the `guess_reader` function that will help guess which
    reader to use from some aspect(s) of the filename.

    (This `guess_reader` function also accepts a `format=`
    keyword that allows the user to explicitly specify that
    the kirk reader should be used.)

    5. Submit a pull request to the github repository for
    this package, so that other folks can use your handy
    new reader too!
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_kirk_fitted_light_curves", "from_kirk_stellar_spectra"]


def from_kirk_fitted_light_curves(self, filepath):
    """
    Populate a Rainbow from a file in James Kirk's fitted
    light curve format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load. This is interpreted to
        be pointing to a directory called `wb_lcs/` that's
        filled with wavelength-binned light curves. It will
        expect there to also be a transmission spectrum
        file called 'transmission_spectrum.data' at the
        same directory level as 'wb_lcs'.

            some-directory/
                wb_lcs/
                transmission_spectrum.dat
    """

    # the directory that contains the light curve
    lc_directory = os.path.abspath(filepath)

    # the expected location of the transmission spectrum
    transmission_spectrum_filepath = lc_directory.replace(
        "wb_lcs", "transmission_spectrum.dat"
    )

    # complain if the spectrum isn't found
    try:
        transmission_spectrum_file = glob.glob(transmission_spectrum_filepath)[0]
    except IndexError:
        raise ValueError(
            f"""
        The `from_kirk_fitted_lightcurves` reader is expecting to find
        {transmission_spectrum_filepath}
        but it's not there 🤨.
        """
        )

    # load the transmission spectrum
    transmission_spectrum = ascii.read(
        transmission_spectrum_file,
        names=[
            "wavelength",
            "wavelength_width",
            "radius_ratio",
            "radius_ratio_uncertainty_upper",
            "radius_ratio_uncertainty_lower",
        ],
    )

    # populate a 1D array of wavelengths (with astropy units of length)
    for k in transmission_spectrum.colnames:
        if "wavelength" in k:
            unit = u.micron
        else:
            unit = 1
        self.wavelike[k] = transmission_spectrum[k].data * unit
    self.wavelike["wavelength_lower"] = (
        self.wavelike["wavelength"] - self.wavelike["wavelength_width"] / 2
    )
    self.wavelike["wavelength_upper"] = (
        self.wavelike["wavelength"] + self.wavelike["wavelength_width"] / 2
    )

    # the filenames for the individual light curve files
    lc_filenames = sorted(glob.glob(os.path.join(lc_directory, "model_tab_wb*.dat")))
    assert len(lc_filenames) == self.nwave

    # loop through the files and load them into arrays (faster than tables)
    lcs = []
    for f in tqdm(lc_filenames):
        lc = np.loadtxt(f)
        lcs.append(lc)
    names = [
        "time",
        "flux",
        "uncertainty",
        "rescaled_uncertainty",
        "model",
        "planet_model",
        "systematics_model",
        "residuals",
    ]

    # extract the time and the fluxlike quantities
    for i, k in enumerate(names):
        if k == "time":
            self.timelike["time"] = lcs[0][:, 0] * u.day
            for lc in lcs:
                assert np.shape(lc)[0] == self.ntime
        else:
            self.fluxlike[k] = np.array([lc[:, i] for lc in lcs])


def from_kirk_stellar_spectra(self, filepath):
    """
    Populate a Rainbow from a file in James Kirk's extracted
    stellar spectra format. This expects a file path pointing
    to a wb_lcs`

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file to load. This is interpreted to
        be pointing to a file called '*flux_resampled_*.pickle',
        and expects there to be at the same directory level two
        other files called '*error_resampled_*.pickle' and
        'wavelengths.pickle'.

            some-directory/
                star1_flux_resampled_unsmoothed.pickle
                star1_erro_resampled_unsmoothed.pickle
                wavelengths.pickle
    """

    # load the flux pickle
    flux_file = filepath
    self.fluxlike["flux"] = pickle.load(open(flux_file, "rb")).T

    # load the uncertainty pickle
    uncertainty_file = filepath.replace("_flux_resampled_", "_error_resampled_")
    self.fluxlike["uncertainty"] = pickle.load(open(uncertainty_file, "rb")).T

    # load the wavelength pickle
    wavelength_file = flux_file.replace(
        os.path.basename(flux_file), "wavelengths.pickle"
    )
    self.wavelike["wavelength"] = pickle.load(open(wavelength_file, "rb")) * u.micron

    warnings.warn(
        f"""
    We're not sure what the times look like for this dataset.
    Making up some totally imaginary ones!
    """
    )
    self.timelike["time"] = np.arange(self.flux.shape[1]) * u.s
