from ...imports import *
from ...resampling import *

__all__ = ["bin", "bin_in_time", "bin_in_wavelength"]


def bin(
    self,
    dt=None,
    time=None,
    ntimes=None,
    R=None,
    dw=None,
    wavelength=None,
    wavelength_edges=None,
    nwavelengths=None,
):
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
    ntimes : int
        A fixed number of times to bin together.
    R : float
        The spectral resolution for creating a grid
        that is uniform in logarithmic space.
    dw : astropy.units.Quantity
        The d(wavelength) bin size for creating a grid
        that is uniform in linear space.
    wavelength : array of astropy.units.Quantity
        An array of wavelengths for the centers of bins,
        if you just want to give an entirely custom array.
    wavelength_edges : array of astropy.units.Quantity
        An array of wavelengths for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `wavelength_edges[:-1]` to
        `wavelength_edges[1:]`, so the resulting binned
        Rainbow will have `len(wavelength_edges) - 1`
        wavelength bins associated with it.
    nwavelengths : int
        A fixed number of wavelengths to bin together.

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """
    # bin first in time
    binned_in_time = self.bin_in_time(dt=dt, time=time, ntimes=ntimes)

    # then bin in wavelength
    binned = binned_in_time.bin_in_wavelength(
        R=R,
        dw=dw,
        wavelength=wavelength,
        wavelength_edges=wavelength_edges,
        nwavelengths=nwavelengths,
    )

    # return the binned object
    return binned


def bin_in_time(self, dt=None, time=None, ntimes=None):
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
    ntimes : int
        A fixed number of times to bin together.

    The time-setting order of precendence is:
        1) time
        2) dt

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    if ntimes is not None:
        raise RuntimeWarning("🌈 Your binning option isn't implemented yet. Sorry! 😔")

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
    # Technically, we should include uncertainties here too,
    # so that times/wavelengths are weighted more toward
    # inputs with higher flux weights (e.g. smaller variance),
    # but that will make non-uniform grids that will be
    # really hard to deal with.
    new.timelike = {}
    for k in self.timelike:
        bt, bv = bintogrid(x=self.time, y=self.timelike[k], unc=None, **binkw)
        new.timelike[k] = bv
    new.timelike[k] = bt

    # bin the flux-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    # TODO (think about cleverer bintogrid for 2D arrays?)
    new.fluxlike = {}
    ok = self.is_ok()
    for k in self.fluxlike:

        if k == "uncertainty":
            warnings.warn(
                """
            Uncertainties and/or data quality flags might
            not be handled absolutely perfectly yet...
            """
            )

        # loop through wavelengths
        for w in range(new.nwave):
            ok_times = ok[w, :]

            if self.uncertainty is None:
                bt, bv = bintogrid(
                    x=self.time[:],
                    y=self.fluxlike[k][w, :],
                    unc=None,
                    **binkw,
                )
                bu = None
            else:
                bt, bv, bu = bintogrid(
                    x=self.time[:],
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

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scale, after binning
    new._guess_wscale()

    return new


def bin_in_wavelength(
    self, R=None, dw=None, wavelength=None, wavelength_edges=None, nwavelengths=None
):
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
    wavelength_edges : array of astropy.units.Quantity
        An array of wavelengths for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `wavelength_edges[:-1]` to
        `wavelength_edges[1:]`, so the resulting binned
        Rainbow will have `len(wavelength_edges) - 1`
        wavelength bins associated with it.
    nwavelengths : int
        A fixed number of wavelengths to bin together.

    The wavelength-setting order of precendence is:
        1) wavelength
        2) dw
        3) R

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    if (nwavelengths is not None) or (wavelength_edges is not None):
        raise RuntimeWarning("🌈 Your binning option isn't implemented yet. Sorry! 😔")

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

    # get a fluxlike array of what's OK to include in the bins
    ok = self.is_ok()
    for k in self.fluxlike:
        # self.speak(f" binning '{k}' in wavelength")
        # self.speak(f"  original shape was {np.shape(self.fluxlike[k])}")
        if k == "uncertainty":
            warnings.warn(
                """
            Uncertainties and/or data quality flags might
            not be handled absolutely perfectly yet...
            """
            )

        for t in range(new.ntime):
            ok_wavelengths = ok[:, t]
            if self.uncertainty is None:
                bt, bv = binning_function(
                    x=self.wavelength[:],
                    y=self.fluxlike[k][:, t],
                    unc=None,
                    **binkw,
                )
                bu = None
            else:
                bt, bv, bu = binning_function(
                    x=self.wavelength[:],
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

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scale, after binning
    new._guess_wscale()
    # new.metadata["wscale"] = wscale
    return new
