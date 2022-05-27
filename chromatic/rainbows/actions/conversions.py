from ...imports import *

__all__ = ["to_nparray", "to_df"]


def to_nparray(self, timeformat='d'):
    """ Convert Rainbow object to 1D and 2D numpy arrays
    Parameters
    ----------
        self : Rainbow object
            chromatic Rainbow object to convert into array format
        timeformat : str
            (optional, default='d')
            The time format to use (seconds, minutes, hours, days etc.)
    Returns
    ----------
        rflux : np.array
            flux (MJy/sr)         [n_wavelengths x n_integrations]
        rfluxu : np.array
            flux uncertainty      [n_wavelengths x n_integrations]
        rtime : np.array
            time (BJD_TDB, hours) [n_integrations]
        rwavel : np.array
            wavelength (microns)  [n_wavelengths]
    """
    secondformat = ['second', 'seconds', 'sec', 's']
    minuteformat = ['minute', 'minutes', 'min', 'm']
    hourformat = ['hour', 'hours', 'h']
    dayformat = ['day', 'days', 'd']
    yearformat = ['year', 'years', 'y']

    rflux  = np.array(self.fluxlike['flux'])  # flux                         : [n_wavelengths x n_integrations]
    rfluxu = np.array(self.fluxlike['uncertainty'])  # uncertainty           : [n_wavelengths x n_integrations]
    rtime  = np.array(self.timelike['time'])  # time (BJD_TDB, hours)        : [n_integrations]
    rwavel  = np.array(self.wavelike['wavelength'])  # wavelength (microns)  : [n_wavelengths]

    # change the time array into the requested format (e.g. seconds, minutes, days etc.)
    if timeformat in secondformat:
        rtime = rtime * 3600
    elif timeformat in minuteformat:
        rtime = rtime * 60
    elif timeformat in hourformat:
        pass
    elif timeformat in dayformat:
        # hours is the default time setting
        rtime = rtime / 24.
    elif timeformat in yearformat:
        rtime = rtime / (24 * 365.)
    else:
        warnings.warn("Unrecognised Time Format!")
        return

    return rflux, rfluxu, rtime, rwavel


def to_df(self, timeformat='d'):
    """ Convert Rainbow object to pandas dataframe
    Parameters
    ----------
        self : Rainbow object
            chromatic Rainbow object to convert into pandas df format
        timeformat : str
            (optional, default='d')
            The time format to use (seconds, minutes, hours, days etc.)
    Returns
    ----------
        pd.DataFrame
    """
    # extract 1D/2D formats for the flux, uncertainty, time and wavelengths
    rflux, rfluxu, rtime, rwavel = self.to_nparray(timeformat)

    # put all arrays onto same dimensions
    x, y = np.meshgrid(rtime, rwavel)
    rainbow_dict = {f"Time ({timeformat})": x.ravel(), "Wavelength (microns)": y.ravel(), "Flux": rflux.ravel(),
                    "Flux Uncertainty": rfluxu.ravel()}

    # convert to pandas dataframe
    df = pd.DataFrame(rainbow_dict)
    return df
