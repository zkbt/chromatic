from ...imports import *

__all__ = ["to_nparray", "to_df"]


def to_nparray(self, t_unit="d", w_unit="micron"):
    """Convert Rainbow object to 1D and 2D numpy arrays
    Parameters
    ----------
        self : Rainbow object
            chromatic Rainbow object to convert into array format
        t_unit : str
            (optional, default='d')
            The time units to use (seconds, minutes, hours, days etc.)
        w_unit : str
            (optional, default='micron')
            The wavelength units to use
    Returns
    ----------
        rflux : array
            flux                   [n_wavelengths x n_integrations]
        rfluxu : array
            flux uncertainty       [n_wavelengths x n_integrations]
        rtime : array
            time (t_unit)          [n_integrations]
        rwavel : array
            wavelength (w_unit)    [n_wavelengths]
    """

    rflux = np.array(
        self.fluxlike["flux"]
    )  # flux                       : [n_wavelengths x n_integrations]
    rfluxu = np.array(
        self.fluxlike["uncertainty"]
    )  # uncertainty        : [n_wavelengths x n_integrations]
    rtime = self.timelike[
        "time"
    ]  # time (BJD_TDB, hours)                : [n_integrations]
    rwavel = self.wavelike[
        "wavelength"
    ]  # wavelength (microns)          : [n_wavelengths]

    try:
        # nice bit of code copied from .imshow(), thanks Zach!
        # Change time into the units requested by the user
        t_unit = u.Unit(t_unit)
    except:
        cheerfully_suggest("Unrecognised Time Format! Returning day by default")
        t_unit = u.Unit("d")
    rtime = np.array(rtime.to(t_unit).value)

    try:
        # Change wavelength into the units requested by the user
        w_unit = u.Unit(w_unit)
    except:
        cheerfully_suggest(
            "Unrecognised Wavelength Format! Returning micron by default"
        )
        w_unit = u.Unit("micron")
    rwavel = np.array(rwavel.to(w_unit).value)

    return rflux, rfluxu, rtime, rwavel


def to_df(self, t_unit="d", w_unit="micron"):
    """Convert Rainbow object to pandas dataframe
    Parameters
    ----------
        self : Rainbow object
            chromatic Rainbow object to convert into pandas df format
        t_unit : str
            (optional, default='d')
            The time units to use (seconds, minutes, hours, days etc.)
        w_unit : str
            (optional, default='micron')
            The wavelength units to use
    Returns
    ----------
        pd.DataFrame
            The rainbow object flattened and converted into a pandas dataframe
    """
    # extract 1D/2D formats for the flux, uncertainty, time and wavelengths
    rflux, rfluxu, rtime, rwavel = self.to_nparray(t_unit=t_unit, w_unit=w_unit)

    # put all arrays onto same dimensions
    x, y = np.meshgrid(rtime, rwavel)
    rainbow_dict = {
        f"Time ({t_unit})": x.ravel(),
        f"Wavelength ({w_unit})": y.ravel(),
        "Flux": rflux.ravel(),
        "Flux Uncertainty": rfluxu.ravel(),
    }

    # convert to pandas dataframe
    df = pd.DataFrame(rainbow_dict)
    return df
