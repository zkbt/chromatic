from ...imports import *
from ...resampling import *

__all__ = ["bin", "bin_in_time", "bin_in_wavelength", "get_integrated_lightcurve"]


def bin(
    self,
    dt=None,
    time=None,
    time_edges=None,
    ntimes=None,
    R=None,
    dw=None,
    wavelength=None,
    wavelength_edges=None,
    nwavelengths=None,
    starting_wavelengths="1D",
    trim=True,
):
    """
    Bin in wavelength and/or time.

    By default, time binning happens before wavelength binning.
    To control the order, use separate calls to `.bin()`.

    The time-setting order of precendence is
    [`time_edges`, `time`, `dt`, `ntimes`]
    The first will be used, and others will be ignored.

    The wavelength-setting order of precendence is
    [`wavelength_edges`, `wavelength`, `dw`, `R`, `nwavelengths`]
    The first will be used, and others will be ignored.


    Parameters
    ----------
    dt : astropy.units.Quantity
        The d(time) bin size for creating a grid
        that is uniform in linear space.
    time : array of astropy.units.Quantity
        An array of times, if you just want to give
        it an entirely custom array.
        The widths of the bins will be guessed from the centers
        (well, if the spacing is uniform constant; pretty well
        but not perfectly otherwise).
    time_edges : array of astropy.units.Quantity
        An array of times for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `time_edges[:-1]` to
        `time_edges[1:]`, so the resulting binned
        Rainbow will have `len(time_edges) - 1`
        time bins associated with it.
    ntimes : int
        A fixed number of time to bin together.
        Binning will start from the 0th element of the
        starting times; if you want to start from
        a different index, trim before binning.
    R : float
        The spectral resolution for creating a grid
        that is uniform in logarithmic space.
    dw : astropy.units.Quantity
        The d(wavelength) bin size for creating a grid
        that is uniform in linear space.
    wavelength : array of astropy.units.Quantity
        An array of wavelengths for the centers of bins,
        if you just want to give an entirely custom array.
        The widths of the bins will be guessed from the centers
        (well, if the spacing is uniform constant; pretty well
        but not perfectly otherwise).
    wavelength_edges : array of astropy.units.Quantity
        An array of wavelengths for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `wavelength_edges[:-1]` to
        `wavelength_edges[1:]`, so the resulting binned
        Rainbow will have `len(wavelength_edges) - 1`
        wavelength bins associated with it.
    nwavelengths : int
        A fixed number of wavelengths to bin together.
        Binning will start from the 0th element of the
        starting wavelengths; if you want to start from
        a different index, trim before binning.
    trim : bool
        Should any wavelengths or columns that end up
        as entirely nan be trimmed out of the result?
        (default = True)

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    # bin first in time
    binned_in_time = self.bin_in_time(
        dt=dt, time=time, time_edges=time_edges, ntimes=ntimes, trim=trim
    )

    # then bin in wavelength
    binned = binned_in_time.bin_in_wavelength(
        R=R,
        dw=dw,
        wavelength=wavelength,
        wavelength_edges=wavelength_edges,
        nwavelengths=nwavelengths,
        starting_wavelengths=starting_wavelengths,
        trim=trim,
    )

    # return the binned object
    return binned


def bin_in_time(self, dt=None, time=None, time_edges=None, ntimes=None, trim=True):
    """
    Bin in time.

    The time-setting order of precendence is
    [`time_edges`, `time`, `dt`, `ntimes`]
    The first will be used, and others will be ignored.


    Parameters
    ----------
    dt : astropy.units.Quantity
        The d(time) bin size for creating a grid
        that is uniform in linear space.
    time : array of astropy.units.Quantity
        An array of times, if you just want to give
        it an entirely custom array.
        The widths of the bins will be guessed from the centers
        (well, if the spacing is uniform constant; pretty well
        but not perfectly otherwise).
    time_edges : array of astropy.units.Quantity
        An array of times for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `time_edges[:-1]` to
        `time_edges[1:]`, so the resulting binned
        Rainbow will have `len(time_edges) - 1`
        time bins associated with it.
    ntimes : int
        A fixed number of time to bin together.
        Binning will start from the 0th element of the
        starting times; if you want to start from
        a different index, trim before binning.
    trim : bool
        Should any wavelengths or columns that end up
        as entirely nan be trimmed out of the result?
        (default = True)

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("bin_in_time", locals())

    # if no bin information is provided, don't bin
    if np.all([x is None for x in [dt, time, time_edges, ntimes]]):
        return self

    # set up binning parameters
    binkw = dict(weighting="inversevariance", drop_nans=False)

    # [`time_edges`, `time`, `dt`, `ntimes`]
    if time_edges is not None:
        binkw["newx_edges"] = time_edges
    elif time is not None:
        binkw["newx"] = time
    elif dt is not None:
        binkw["dx"] = dt
    elif ntimes is not None:
        binkw["nx"] = ntimes

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
        binned = bintogrid(x=self.time, y=self.timelike[k], unc=None, **binkw)
        new.timelike[k] = binned["y"]
    new.timelike["time"] = binned["x"]
    new.timelike["time_lower"] = binned["x_edge_lower"]
    new.timelike["time_upper"] = binned["x_edge_upper"]

    # bin the flux-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    # TODO (think about cleverer bintogrid for 2D arrays?)
    new.fluxlike = {}
    ok = self.is_ok()
    for k in self.fluxlike:

        '''
        if k == "uncertainty":
            warnings.warn(
                """
            Uncertainties and/or data quality flags might
            not be handled absolutely perfectly yet...
            """
            )'''

        # loop through wavelengths
        for w in range(new.nwave):

            # FIXME - not being used yet?
            ok_times = ok[w, :]

            # set what uncertainties should be used for binning
            if self.uncertainty is None:
                uncertainty_for_binning = None
            else:
                uncertainty_for_binning = self.uncertainty[w, :]

            # bin the quantities for this wavelength
            binned = bintogrid(
                x=self.time[:],
                y=self.fluxlike[k][w, :],
                unc=uncertainty_for_binning,
                **binkw,
            )

            # if necessary, create a new fluxlike array
            if k not in new.fluxlike:
                new_shape = (new.nwave, new.ntime)
                new.fluxlike[k] = np.zeros(new_shape)
                # FIXME make this more robust to units

            # store the binned array in the appropriate place
            if k == "uncertainty":
                # uncertainties are usually standard error on the mean
                new.fluxlike[k][w, :] = binned["uncertainty"]
            else:
                # note: all quantities are weighted the same as flux (probably inversevariance)
                new.fluxlike[k][w, :] = binned["y"]

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scale, after binning
    new._guess_wscale()

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow (with trimming if necessary)
    if trim:
        return new.trim_nan_times()
    else:
        return new


def bin_in_wavelength(
    self,
    R=None,
    dw=None,
    wavelength=None,
    wavelength_edges=None,
    nwavelengths=None,
    starting_wavelengths="1D",
    trim=True,
):
    """
    Bin in wavelength.

    The wavelength-setting order of precendence is
    [`wavelength_edges`, `wavelength`, `dw`, `R`, `nwavelengths`]
    The first will be used, and others will be ignored.

    Parameters
    ----------
    R : float
        The spectral resolution for creating a grid
        that is uniform in logarithmic space.
    dw : astropy.units.Quantity
        The d(wavelength) bin size for creating a grid
        that is uniform in linear space.
    wavelength : array of astropy.units.Quantity
        An array of wavelength centers, if you just want to give
        it an entirely custom array. The widths of the bins
        will be guessed from the centers. It will do a good
        job if the widths are constant, but don't 100% trust
        it otherwise.
    wavelength_edges : array of astropy.units.Quantity
        An array of wavelengths for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `wavelength_edges[:-1]` to
        `wavelength_edges[1:]`, so the resulting binned
        Rainbow will have `len(wavelength_edges) - 1`
        wavelength bins associated with it.
    nwavelengths : int
        A fixed number of wavelengths to bin together.
        Binning will start from the 0th element of the
        starting wavelengths; if you want to start from
        a different index, trim before binning.
    trim : bool
        Should any wavelengths or columns that end up
        as entirely nan be trimmed out of the result?
        (default = True)

    Returns
    -------
    binned : Rainbow
        The binned Rainbow.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("bin_in_wavelength", locals())

    # if no bin information is provided, don't bin
    if (
        (wavelength is None)
        and (wavelength_edges is None)
        and (nwavelengths is None)
        and (dw is None)
        and (R is None)
    ):
        return self

    # set up binning parameters
    binkw = dict(weighting="inversevariance", drop_nans=False)

    # [`wavelength_edges`, `wavelength`, `dw`, `R`, `nwavelengths`]
    if wavelength_edges is not None:
        binning_function = bintogrid
        binkw["newx_edges"] = wavelength_edges
    elif wavelength is not None:
        binning_function = bintogrid
        binkw["newx"] = wavelength
    elif dw is not None:
        binning_function = bintogrid
        binkw["dx"] = dw
    elif R is not None:
        binning_function = bintoR
        binkw["R"] = R
    elif nwavelengths is not None:
        binning_function = bintogrid
        binkw["nx"] = nwavelengths

    # create a new, empty Rainbow
    new = self._create_copy()

    # populate the time information
    new.timelike = {**self.timelike}

    # bin the time-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    new.wavelike = {}
    for k in self.wavelike:
        binned = binning_function(
            x=self.wavelength, y=self.wavelike[k], unc=None, **binkw
        )
        new.wavelike[k] = binned["y"]
    new.wavelike["wavelength"] = binned["x"]
    new.wavelike["wavelength_lower"] = binned["x_edge_lower"]
    new.wavelike["wavelength_upper"] = binned["x_edge_upper"]

    # bin the flux-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    # TODO (think about cleverer bintogrid for 2D arrays)
    new.fluxlike = {}

    # get a fluxlike array of what's OK to include in the bins
    ok = self.is_ok()
    for k in self.fluxlike:

        '''
        if k == "uncertainty":
            warnings.warn(
                """
            Uncertainties and/or data quality flags might
            not be handled absolutely perfectly yet...
            """
            )'''

        for t in range(new.ntime):

            # FIXME - not being used yet?
            ok_wavelengths = ok[:, t]

            # set what uncertainties should be used for binning
            if self.uncertainty is None:
                uncertainty_for_binning = None
            else:
                uncertainty_for_binning = self.uncertainty[:, t]

            if starting_wavelengths.upper() == "1D":
                w = self.wavelength[:]
            elif starting_wavelengths.upper() == "2D":
                w = self.fluxlike["wavelength"][:, t]
            # bin the quantities for this time
            binned = binning_function(
                x=w,
                y=self.fluxlike[k][:, t],
                unc=uncertainty_for_binning,
                **binkw,
            )

            # if necessary, create a new fluxlike array
            if k not in new.fluxlike:
                new_shape = (new.nwave, new.ntime)
                new.fluxlike[k] = np.zeros(new_shape)
                # FIXME make this more robust to units

            # store the binned array in the appropriate place
            if k == "uncertainty":
                # uncertainties are usually standard error on the mean
                new.fluxlike[k][:, t] = binned["uncertainty"]
            else:
                # note: all quantities are weighted the same as flux (probably inversevariance)
                new.fluxlike[k][:, t] = binned["y"]

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scale, after binning
    new._guess_wscale()

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow (with trimming if necessary)
    if trim:
        return new.trim_nan_wavelengths()
    else:
        return new


def get_integrated_lightcurve(self):
    """
    Produce a wavelength-integrated light curve.

    Returns
    -------
    lc : Rainbow
        A Rainbow object with just one wavelength.
    """
    return self.bin(nwavelengths=self.nwave)
