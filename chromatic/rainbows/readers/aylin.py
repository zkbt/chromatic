# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_aylin"]


def from_aylin(self, filepaths):
    """
    Populate a Rainbow from a file in the aylin format.

    Parameters
    ----------

    filepaths : str, list
        A glob-friendly string pointing to all the spectrum files
        for a particular target, such as 'Halpha/spec1d*Original_Spectrum.csv'.
            ...or...
        alteternatively, a list of 'Original_Spectrum.csv' files.

        In both cases, the code will assume there are parallel 'Shifted_Wave.csv'
        for all 'Original_Spectrum' files, and use them to load wavelengths
        with the H-alpha line aligned.
    """

    if isinstance(filepaths, str):
        files = glob.glob(filepaths)
    elif isinstance(filepaths, list):
        files = filepaths
    else:
        raise ValueError(
            """
        The input to the 'read_aylin' reader should be a string
        with a wildcard in it like 'spec1d*Original_Spectrum.csv',
        or a list of files like that.
        """
        )

    # load all the extracted spectrum files
    tables = {}
    for f in files:
        tables[f] = ascii.read(f)

    # figure out the number of times and maximum number of wavelengths
    N_times = len(files)
    N_wavelengths = np.max([len(t) for t in tables.values()])

    # create empty 2D arrays
    unshifted_wavelengths = np.zeros((N_wavelengths, N_times)) * np.nan
    shifted_wavelengths = np.zeros((N_wavelengths, N_times)) * np.nan
    fluxes = np.zeros((N_wavelengths, N_times)) * np.nan
    flux_uncertainties = np.zeros((N_wavelengths, N_times)) * np.nan

    # create empty list for times
    times = []

    # loop through spectra
    for i, f in enumerate(tables):

        # get this spectrum
        t = tables[f]

        # use only finite wavelengths
        N_ok = len(t)
        unshifted_wavelengths[:N_ok, i] = t["wavelength"]
        fluxes[:N_ok, i] = t["flux"]
        flux_uncertainties[:N_ok, i] = t["flux_err"]

        # load the shifted wavelengths
        shifted_wavelength_filename = f.replace(
            "Original_Spectrum.csv", "Shifted_Wave.csv"
        )
        wavelength_table = ascii.read(
            shifted_wavelength_filename,
            delimiter=",",
            data_start=1,
            names=["index", "wavelength"],
        )
        shifted_wavelengths[:N_ok, i] = wavelength_table["wavelength"]

        # figure out the time from the filename
        s = f.split("Original_Spectrum")[0].split("_")[-1]
        time = (
            Time(
                {
                    "year": int(s[0:4]),
                    "month": int(s[4:6]),
                    "day": int(s[6:8]),
                    "hour": int(s[9:11]),
                    "minute": int(s[11:13]),
                    "second": float(s[13:]),
                },
                scale="utc",
            ).jd
            * u.day
        )
        times.append(time)

    self.fluxlike["flux"] = fluxes
    self.fluxlike["uncertainty"] = flux_uncertainties
    self.fluxlike["unshifted_wavelength_2d"] = unshifted_wavelengths * u.Angstrom
    self.fluxlike["wavelength_2d"] = shifted_wavelengths * u.Angstrom

    self.wavelike["wavelength"] = np.mean(self.fluxlike["wavelength_2d"], axis=1)
    self.timelike["time"] = u.Quantity(times)
    # self.timelike["filename"] = np.array(files)

    cheerfully_suggest(
        """
    These spectra may have slightly different wavelength grids for different time-points.
    These per-time wavelength arrays are stored in rainbow.fluxlike['wavelength_2d'].
    To align the fluxes onto a shared wavelength grid, we strongly suggest running
        `aligned = rainbow.align_wavelengths()`
    which will shift the flux values so that everything is on one shared wavelength grid,
    at the cost of introducing slight correlations between adjacent data points.
    """
    )
