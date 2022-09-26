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
    for f in tqdm(lc_filenames, leave=False):
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
                star1_error_resampled_unsmoothed.pickle
                wvl_solution.pickle
                BJD_TDB_time.pickle
    """

    # load the flux pickle
    flux_file = filepath
    self.fluxlike["flux"] = pickle.load(open(flux_file, "rb")).T

    # load the uncertainty pickle
    uncertainty_file = filepath.replace("_flux_resampled", "_error_resampled")
    assert uncertainty_file != flux_file
    self.fluxlike["uncertainty"] = pickle.load(open(uncertainty_file, "rb")).T

    # load the wavelength pickle
    wavelength_file = flux_file.replace(
        os.path.basename(flux_file), "wvl_solution.pickle"
    )
    self.wavelike["wavelength"] = pickle.load(open(wavelength_file, "rb")) * u.micron

    # load the wavelength pickle
    time_file = flux_file.replace(os.path.basename(flux_file), "BJD_TDB_time.pickle")
    bjd_mjd = pickle.load(open(time_file, "rb"))
    astropy_times = Time(bjd_mjd, format="mjd", scale="tdb")
    self.set_times_from_astropy(astropy_times, is_barycentric=True)
