from ..imports import *
from ..resampling import *


def normalize(self, wavelength=True, time=False):
    """
    Normalize by dividing through by the median spectrum and/or lightcurve.

    Parameters
    ----------
    wavelength : bool
        Should we divide by the median spectrum?

    time : bool
        Should we divide by the median light curve?

    Returns
    -------
    normalized : MultiRainbow
        The normalized MultiRainbow.
    """

    # TODO, think about more careful treatment of uncertainties + good/bad data
    new = self._create_copy()

    # (ignore nan warnings)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if time:
            new = new / np.nanmedian(self.flux, axis=self.waveaxis)

        if wavelength:
            new = new / np.nanmedian(self.flux, axis=self.timeaxis)

    return new


def bin(self, dt=None, time=None, R=None, dw=None, wavelength=None):
    """
    Bin the rainbow in wavelength and/or time.

    The time-setting order of precendence is
    [`time`, `dt`], meaning that if `time` is set,
    any values given for `dt` will be ignored.
    The wavelength-setting order of precendence is
    [`wavelength`, `dw`, `R`], meaning that if `wavelength`
    is set any values of `dw` or `R` will be ignored, and
    if `dw` is set any value of `R` will be ignored.

    Parameters
    ----------
    dt : astropy.units.Quantity
        The d(time) bin size for creating a grid
        that is uniform in linear space.
    time : array of astropy.units.Quantity
        An array of times, if you just want to give
        it an entirely custom array.
    R : float
        The spectral resolution for creating a grid
        that is uniform in logarithmic space.
    dw : astropy.units.Quantity
        The d(wavelength) bin size for creating a grid
        that is uniform in linear space.
    wavelength : array of astropy.units.Quantity
        An array of wavelengths, if you just want to give
        it an entirely custom array.

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    # self.speak(f'binning')

    # bin first in time
    binned_in_time = self.bin_in_time(dt=dt, time=time)

    # then bin in wavelength
    binned = binned_in_time.bin_in_wavelength(R=R, dw=dw, wavelength=wavelength)

    # return the binned object
    return binned


def bin_in_time(self, dt=None, time=None):
    """
    Bin the rainbow in time.

    Parameters
    ----------
    dt : astropy.units.Quantity
        The d(time) bin size for creating a grid
        that is uniform in linear space.
    time : array of astropy.units.Quantity
        An array of times, if you just want to give
        it an entirely custom array.

    The time-setting order of precendence is:
        1) time
        2) dt

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    # if no bin information is provided, don't bin
    if (time is None) and (dt is None):
        return self

    # set up binning parameters
    binkw = dict(weighting="inversevariance", drop_nans=False)
    # TODO: make sure drop_nans doesn't skip rows or columns unexpectedly
    if time is not None:
        binkw["newx"] = time
        # self.speak(f'binning to time={time}')
    elif dt is not None:
        binkw["dx"] = dt
        # self.speak(f'binning to dt={dt}')

    # create a new, empty Rainbow
    new = self._create_copy()

    # populate the wavelength information
    new.wavelike = {**self.wavelike}
    new.metadata["wscale"] = self.wscale

    # bin the time-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    new.timelike = {}
    for k in self.timelike:
        bt, bv = bintogrid(x=self.time, y=self.timelike[k], unc=None, **binkw)
        new.timelike[k] = bv
    new.timelike[k] = bt

    # bin the flux-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    # TODO (think about cleverer bintogrid for 2D arrays)
    new.fluxlike = {}
    for k in self.fluxlike:
        # self.speak(f" binning '{k}' in time")
        # self.speak(f"  original shape was {np.shape(self.fluxlike[k])}")

        if k == "uncertainty":
            warnings.warn(
                """
            Uncertainties might not be handled well yet...
            """
            )
            continue

        for w in range(new.nwave):
            if self.uncertainty is None:
                bt, bv = bintogrid(
                    x=self.time, y=self.fluxlike[k][w, :], unc=None, **binkw
                )
                bu = None
            else:
                bt, bv, bu = bintogrid(
                    x=self.time,
                    y=self.fluxlike[k][w, :],
                    unc=self.uncertainty[w, :],
                    **binkw,
                )

            if k not in new.fluxlike:
                new_shape = (new.nwave, new.ntime)
                new.fluxlike[k] = np.zeros(new_shape)
                # TODO make this more robust to units

            if k == "uncertainty":
                new.fluxlike[k][w, :] = bu
            else:
                new.fluxlike[k][w, :] = bv
        # self.speak(f"  new shape is {np.shape(new.fluxlike[k])}")

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scale, after binning
    new._guess_wscale()

    return new


def bin_in_wavelength(self, R=None, dw=None, wavelength=None):
    """
    Bin the rainbow in wavelength.

    Parameters
    ----------
    R : float
        The spectral resolution for creating a grid
        that is uniform in logarithmic space.
    dw : astropy.units.Quantity
        The d(wavelength) bin size for creating a grid
        that is uniform in linear space.
    wavelength : array of astropy.units.Quantity
        An array of wavelengths, if you just want to give
        it an entirely custom array.

    The wavelength-setting order of precendence is:
        1) wavelength
        2) dw
        3) R

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    # if no bin information is provided, don't bin
    if (wavelength is None) and (dw is None) and (R is None):
        return self

    # set up binning parameters
    binkw = dict(weighting="inversevariance", drop_nans=False)
    if wavelength is not None:
        binning_function = bintogrid
        binkw["newx"] = wavelength
        # self.speak(f'binning to wavelength={wavelength}')
    elif dw is not None:
        binning_function = bintogrid
        binkw["dx"] = dw
        # self.speak(f'binning to dw={dw}')
    elif R is not None:
        binning_function = bintoR
        binkw["R"] = R
        # self.speak(f'binning to R={R}')

    # create a new, empty Rainbow
    new = self._create_copy()

    # populate the time information
    new.timelike = {**self.timelike}

    # bin the time-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    new.wavelike = {}
    for k in self.wavelike:
        bt, bv = binning_function(
            x=self.wavelength, y=self.wavelike[k], unc=None, **binkw
        )
        new.wavelike[k] = bv
    new.wavelike[k] = bt

    # bin the flux-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    # TODO (think about cleverer bintogrid for 2D arrays)
    new.fluxlike = {}
    for k in self.fluxlike:
        # self.speak(f" binning '{k}' in wavelength")
        # self.speak(f"  original shape was {np.shape(self.fluxlike[k])}")
        if k == "uncertainty":
            warnings.warn(
                """
            Uncertainties might not be handled well yet...
            """
            )
            continue

        for t in range(new.ntime):
            if self.uncertainty is None:
                bt, bv = binning_function(
                    x=self.wavelength, y=self.fluxlike[k][:, t], unc=None, **binkw
                )
                bu = None
            else:
                bt, bv, bu = binning_function(
                    x=self.wavelength,
                    y=self.fluxlike[k][:, t],
                    unc=self.uncertainty[:, t],
                    **binkw,
                )

            if k not in new.fluxlike:
                new_shape = (new.nwave, new.ntime)
                new.fluxlike[k] = np.zeros(new_shape)
                # TODO make this more robust to units

            if k == "uncertainty":
                new.fluxlike[k][:, t] = bu
            else:
                new.fluxlike[k][:, t] = bv
        # self.speak(f"  new shape is {np.shape(new.fluxlike[k])}")

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scale, after binning
    new._guess_wscale()
    # new.metadata["wscale"] = wscale
    return new
