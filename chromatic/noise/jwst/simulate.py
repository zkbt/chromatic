from ...imports import *
from ...rainbows import SimulatedRainbow


def make_simulated_jwst_rainbow(table, dt=None, tlim=None, which_snr_to_use=None):
    """
    Make a simulated Rainbow object with noise and cadence set by PandExo.

    Parameters
    ----------
    table : Table
        The 1D tabular output from either PandExo or ETC.
    dt : Quantity
        The time cadence, with units of time.
        (default of None will get integration time from PandExo)
    tlim : Quantity
        The time span of the observation [tmin, tmax], with units of time.
        (default of None will set to two transit durations)
    which_snr_to_use : str
        Which choice of `snr_per_integration_...` should be used?
    """

    # set the wavelengths
    wavelength = table["wavelength"] * u.micron

    # set the cadence (or use from PandExo)
    if dt is None:
        dt = table.meta["time_per_integration"] * u.s

    # set the duration of the observation
    if tlim is None:
        transit_duration = table.meta["transit_duration"] * u.hour
        tlim = transit_duration * [-1, 1]

    # set the S/N per integration per pixel
    if which_snr_to_use is None:
        snr = np.inf
        for k in table.colnames:
            if "snr_per_integration" in k:
                this_snr = np.median(table[k])
                if this_snr < snr:
                    snr = this_snr
                    which_snr_to_use = k
        print(
            f"Using {k} because its median S/N ({np.median(snr):.0f}) is most cautious."
        )

    # create the simulated rainbow
    s = SimulatedRainbow(tlim=tlim, dt=dt, wavelength=wavelength).inject_noise(
        signal_to_noise=table[which_snr_to_use]
    )
    return s
