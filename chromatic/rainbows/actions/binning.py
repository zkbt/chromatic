from ...imports import *
from ...resampling import *

__all__ = [
    "bin",
    "bin_in_time",
    "bin_in_wavelength",
    "get_average_lightcurve_as_rainbow",
    "get_average_spectrum_as_rainbow",
]


def _warn_about_weird_binning(N, dimension, fraction_that_can_be_bad=0.0):
    """
    A helper to warn about too many bins containing fewer than one
    original data point. This is a warning that should generally
    be raised to someone who is trying to create underpopulated
    bins and hasn't made an explicit decision about what to do
    with them.

    Parameters
    ----------
    N : array
        The effective number of original bins going into each new bin.
    minimum_points_per_bin : float
        The threshold for the number of original bins to not be enough.
    dimension : string
        The name of the dimension to warn about = 'wavelength' or 'time'
    """
    if len(N) <= 2:
        return

    N_not_edges = N[1:-1]
    fraction_that_are_bad = np.sum(N_not_edges < 1) / len(N_not_edges)
    if fraction_that_are_bad > fraction_that_can_be_bad:
        message = f"""
        Of the {len(N_not_edges)} non-edge new {dimension} bins,
        {fraction_that_are_bad:.1%} of them effectively contain fewer
        than one original {dimension}.

        Here are your options:
        1) Rerun you binning with larger {dimension} bin sizes to
        decrease the chances that they will be partially populated.
        2) Rerun your binning but change `minimum_points_per_bin=` to a number.
        This will set a lower limit on the effective number of inputs
        points required for each bin. Bins that don't meet this limit
        will be marked as not `ok`, and if `trim=True` (default)
        these bins will automatically be trimmed away. A threshold
        of 1 means bins should average together one or more input
        data; a threshold of 0 will get rid of this warning but
        allow many bins to come from the same data point, so you
        should expect weird correlations.
        """
        cheerfully_suggest(message)


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
    minimum_acceptable_ok=1,
    minimum_points_per_bin=None,
    trim=True,
):
    """
    Bin in wavelength and/or time.

    Average together some number of adjacent data points,
    in wavelength and/or time. For well-behaved data where
    data points are independent from each other, binning down
    by N data points should decrease the noise per bin by
    approximately 1/sqrt(N), making it easier to see subtle
    signals. To bin data points together, data are combined
    using inverse-variance weighting through interpolation
    of cumulative distributions, in an attempt to make sure
    that flux integrals between limits are maintained.

    Currently, the inverse-variance weighting is most reliable
    only for datasets that have been normalized to be close
    to 1. We still need to do a little work to make sure
    it works well on unnormalized datasets with dramatically
    non-uniform uncertainties.

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
    dt : Quantity
        The d(time) bin size for creating a grid
        that is uniform in linear space.
    time : Quantity
        An array of times, if you just want to give
        it an entirely custom array.
        The widths of the bins will be guessed from the centers
        (well, if the spacing is uniform constant; pretty well
        but not perfectly otherwise).
    time_edges : Quantity
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
    dw : Quantity
        The d(wavelength) bin size for creating a grid
        that is uniform in linear space.
    wavelength : Quantity
        An array of wavelengths for the centers of bins,
        if you just want to give an entirely custom array.
        The widths of the bins will be guessed from the centers
        (well, if the spacing is uniform constant; pretty well
        but not perfectly otherwise).
    wavelength_edges : Quantity
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
    minimum_acceptable_ok : float
        The numbers in the `.ok` attribute express "how OK?" each
        data point is, ranging from 0 (not OK) to 1 (super OK).
        In most cases, `.ok` will be binary, but there may be times
        where it's intermediate (for example, if a bin was created
        from some data that were not OK and some that were).
        The `minimum_acceptable_ok` parameter allows you to specify what
        level of OK-ness for a point to go into the binning.
        Reasonable options may include:
            minimum_acceptable_ok = 1
                  Only data points that are perfectly OK
                  will go into the binning.
            minimum_acceptable_ok = 1e-10
                  All data points that aren't definitely not OK
                  will go into the binning.
            minimum_acceptable_ok = 0
                  All data points will be included in the bin.
    minimum_points_per_bin : float
        If you're creating bins that are smaller than those in
        the original dataset, it's possible to end up with bins
        that effectively contain fewer than one original datapoint
        (in the sense that the contribution of one original datapoint
        might be split across multiple new bins). By default,
        we allow this behavior with `minimum_points_per_bin=0`, but you can
        limit your result to only bins that contain one or more
        original datapoints with `minimum_points_per_bin=1`.
    trim : bool
        Should any wavelengths or columns that end up
        as entirely nan be trimmed out of the result?
        (default = True)

    Returns
    -------
    binned : Rainbow
        The binned `Rainbow`.
    """

    # bin first in time
    binned_in_time = self.bin_in_time(
        dt=dt,
        time=time,
        time_edges=time_edges,
        ntimes=ntimes,
        minimum_acceptable_ok=minimum_acceptable_ok,
        minimum_points_per_bin=minimum_points_per_bin,
        trim=trim,
    )

    # then bin in wavelength
    binned = binned_in_time.bin_in_wavelength(
        R=R,
        dw=dw,
        wavelength=wavelength,
        wavelength_edges=wavelength_edges,
        nwavelengths=nwavelengths,
        minimum_acceptable_ok=minimum_acceptable_ok,
        minimum_points_per_bin=minimum_points_per_bin,
        trim=trim,
    )

    # return the binned object
    return binned


def bin_in_time(
    self,
    dt=None,
    time=None,
    time_edges=None,
    ntimes=None,
    minimum_acceptable_ok=1,
    minimum_points_per_bin=None,
    trim=True,
):
    """
    Bin in time.

    The time-setting order of precendence is
    [`time_edges`, `time`, `dt`, `ntimes`]
    The first will be used, and others will be ignored.


    Parameters
    ----------
    dt : Quantity
        The d(time) bin size for creating a grid
        that is uniform in linear space.
    time : Quantity
        An array of times, if you just want to give
        it an entirely custom array.
        The widths of the bins will be guessed from the centers
        (well, if the spacing is uniform constant; pretty well
        but not perfectly otherwise).
    time_edges : Quantity
        An array of times for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `time_edges[:-1]` to
        `time_edges[1:]`, so the resulting binned
        `Rainbow` will have `len(time_edges) - 1`
        time bins associated with it.
    ntimes : int
        A fixed number of time to bin together.
        Binning will start from the 0th element of the
        starting times; if you want to start from
        a different index, trim before binning.
    minimum_acceptable_ok : float
        The numbers in the `.ok` attribute express "how OK?" each
        data point is, ranging from 0 (not OK) to 1 (super OK).
        In most cases, `.ok` will be binary, but there may be times
        where it's intermediate (for example, if a bin was created
        from some data that were not OK and some that were).
        The `minimum_acceptable_ok` parameter allows you to specify what
        level of OK-ness for a point to go into the binning.
        Reasonable options may include:
            minimum_acceptable_ok = 1
                  Only data points that are perfectly OK
                  will go into the binning.
            minimum_acceptable_ok = 1e-10
                  All data points that aren't definitely not OK
                  will go into the binning.
            minimum_acceptable_ok = 0
                  All data points will be included in the bin.
    minimum_points_per_bin : float
        If you're creating bins that are smaller than those in
        the original dataset, it's possible to end up with bins
        that effectively contain fewer than one original datapoint
        (in the sense that the contribution of one original datapoint
        might be split across multiple new bins). By default,
        we allow this behavior with `minimum_points_per_bin=0`, but you can
        limit your result to only bins that contain one or more
        original datapoints with `minimum_points_per_bin=1`.
    trim : bool
        Should any wavelengths or columns that end up
        as entirely nan be trimmed out of the result?
        (default = True)

    Returns
    -------
    binned : Rainbow
        The binned `Rainbow`.
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
    new.timelike["unbinned_times_per_binned_time"] = binned["N_unbinned/N_binned"]

    # bin the flux-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    # TODO (think about cleverer bintogrid for 2D arrays?)
    new.fluxlike = {}
    ok = self.ok
    # loop through wavelengths
    for w in tqdm(np.arange(new.nwave), leave=False):

        '''
        if k == "uncertainty":
            cheerfully_suggest(
                """
            Uncertainties and/or data quality flags might
            not be handled absolutely perfectly yet...
            """
            )'''

        for k in self.fluxlike:
            # mask out "bad" wavelengths
            time_is_bad = ok[w, :] < minimum_acceptable_ok
            if (self.uncertainty is None) or np.all(self.uncertainty == 0):
                uncertainty_for_binning = np.ones(self.ntime).astype(bool)
            elif k in self._keys_that_get_uncertainty_weighting:
                uncertainty_for_binning = self.uncertainty[w, :] * 1
            else:
                uncertainty_for_binning = np.ones(self.ntime).astype(bool)

            if k != "ok":
                uncertainty_for_binning[time_is_bad] = np.inf

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
                if isinstance(self.fluxlike[k], u.Quantity):
                    new.fluxlike[k] *= self.fluxlike[k].unit

            # store the binned array in the appropriate place
            if k == "uncertainty":
                # uncertainties are usually standard error on the mean
                new.fluxlike[k][w, :] = binned["uncertainty"]
            else:
                # note: all quantities are weighted the same as flux (probably inversevariance)
                new.fluxlike[k][w, :] = binned["y"]

    if (new.nwave == 0) or (new.ntime == 0):
        message = f"""
        You tried to bin {self} to {new}.

        After accounting for `minimum_acceptable_ok > {minimum_acceptable_ok}`,
        all new bins would end up with no usable data points.
        Please (a) make sure your input `Rainbow` has at least
        one wavelength and time, (b) check `.ok` accurately expresses
        which data you think are usable, (c) change the `minimum_acceptable_ok`
        keyword for `.bin` to a smaller value, and/or (d) try larger bins.
        """
        cheerfully_suggest(message)
        raise RuntimeError("No good data to bin! (see above)")

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scales, after binning
    new._guess_wscale()
    new._guess_tscale()

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # deal with bins that are smaller than original
    N = new.timelike["unbinned_times_per_binned_time"]
    if minimum_points_per_bin is None:
        _warn_about_weird_binning(N, "time")
    else:
        ok = new.timelike.get("ok", np.ones(new.ntime, bool))
        new.timelike["ok"] = ok * (N >= minimum_points_per_bin)

    # return the new Rainbow (with trimming if necessary)
    if trim:
        return new.trim_times(minimum_acceptable_ok=minimum_acceptable_ok)
    else:
        return new


def bin_in_wavelength(
    self,
    R=None,
    dw=None,
    wavelength=None,
    wavelength_edges=None,
    nwavelengths=None,
    minimum_acceptable_ok=1,
    minimum_points_per_bin=None,
    trim=True,
    starting_wavelengths="1D",
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
    dw : Quantity
        The d(wavelength) bin size for creating a grid
        that is uniform in linear space.
    wavelength : Quantity
        An array of wavelength centers, if you just want to give
        it an entirely custom array. The widths of the bins
        will be guessed from the centers. It will do a good
        job if the widths are constant, but don't 100% trust
        it otherwise.
    wavelength_edges : Quantity
        An array of wavelengths for the edges of bins,
        if you just want to give an entirely custom array.
        The bins will span `wavelength_edges[:-1]` to
        `wavelength_edges[1:]`, so the resulting binned
        `Rainbow` will have `len(wavelength_edges) - 1`
        wavelength bins associated with it.
    nwavelengths : int
        A fixed number of wavelengths to bin together.
        Binning will start from the 0th element of the
        starting wavelengths; if you want to start from
        a different index, trim before binning.
    minimum_acceptable_ok : float
        The numbers in the `.ok` attribute express "how OK?" each
        data point is, ranging from 0 (not OK) to 1 (super OK).
        In most cases, `.ok` will be binary, but there may be times
        where it's intermediate (for example, if a bin was created
        from some data that were not OK and some that were).
        The `minimum_acceptable_ok` parameter allows you to specify what
        level of OK-ness for a point to go into the binning.
        Reasonable options may include:
            minimum_acceptable_ok = 1
                  Only data points that are perfectly OK
                  will go into the binning.
            minimum_acceptable_ok = 1e-10
                  All data points that aren't definitely not OK
                  will go into the binning.
            minimum_acceptable_ok = 0
                  All data points will be included in the bin.
    minimum_points_per_bin : float
        If you're creating bins that are smaller than those in
        the original dataset, it's possible to end up with bins
        that effectively contain fewer than one original datapoint
        (in the sense that the contribution of one original datapoint
        might be split across multiple new bins). By default,
        we allow this behavior with `minimum_points_per_bin=0`, but you can
        limit your result to only bins that contain one or more
        original datapoints with `minimum_points_per_bin=1`.
    trim : bool
        Should any wavelengths or columns that end up
        as entirely nan be trimmed out of the result?
        (default = True)
    starting_wavelengths : str
        What wavelengths should be used as the starting
        value from which we will be binning? Options include:
        '1D' = (default) the shared 1D wavelengths for all times
               stored in `.wavelike['wavelength']`
        '2D' = (used only by `align_wavelengths`) the per-time 2D array
               stored in `.fluxlike['wavelength']`
        [Most users probably don't need to change this from default.]

    Returns
    -------
    binned : Rainbow
        The binned `Rainbow`.
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

    if (
        (self._is_probably_normalized() == False)
        and (self.uncertainty is not None)
        and np.any(self.uncertainty != 0)
    ):

        cheerfully_suggest(
            f"""
        It looks like you're trying to bin in wavelength for a
        `Rainbow` object that might not be normalized. In the
        current version of `chromatic`, binning before normalizing
        might give inaccurate results if the typical uncertainty
        varies strongly with wavelength.

        Please consider normalizing first, for example with
        `rainbow.normalize().bin(...)`
        so that all uncertainties will effectively be relative,
        and the inverse variance weighting used for binning
        wavelengths together will give more reasonable answers.

        If you really need to bin before normalizing, please submit
        an Issue at github.com/zkbt/chromatic/, and we'll try to
        prioritize implementing a statistically sound solution as
        soon as possible!
        """
        )

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
            x=self.wavelike["wavelength"], y=self.wavelike[k], unc=None, **binkw
        )
        new.wavelike[k] = binned["y"]
    new.wavelike["wavelength"] = binned["x"]
    new.wavelike["wavelength_lower"] = binned["x_edge_lower"]
    new.wavelike["wavelength_upper"] = binned["x_edge_upper"]
    new.wavelike["unbinned_wavelengths_per_binned_wavelength"] = binned[
        "N_unbinned/N_binned"
    ]

    # bin the flux-like variables
    # TODO (add more careful treatment of uncertainty + DQ)
    # TODO (think about cleverer bintogrid for 2D arrays)
    new.fluxlike = {}

    # get a fluxlike array of what's OK to include in the bins
    ok = self.ok
    for t in tqdm(np.arange(new.ntime), leave=False):

        for k in self.fluxlike:

            # mask out "bad" wavelengths
            wavelength_is_bad = ok[:, t] < minimum_acceptable_ok

            if (self.uncertainty is None) or np.all(self.uncertainty == 0):
                uncertainty_for_binning = np.ones(self.nwave).astype(bool)
            elif k in self._keys_that_get_uncertainty_weighting:
                uncertainty_for_binning = self.uncertainty[:, t] * 1
            else:
                uncertainty_for_binning = np.ones(self.nwave).astype(bool)
            if k != "ok":
                uncertainty_for_binning[wavelength_is_bad] = np.inf

            if starting_wavelengths.upper() == "1D":
                w = self.wavelike["wavelength"][:]
            elif starting_wavelengths.upper() == "2D":
                w = self.fluxlike["wavelength_2d"][:, t]
            # bin the quantities for this time
            binned = binning_function(
                x=w,
                y=self.fluxlike[k][:, t] * 1,
                unc=uncertainty_for_binning,
                **binkw,
            )

            # if necessary, create a new fluxlike array
            if k not in new.fluxlike:
                new_shape = (new.nwave, new.ntime)
                new.fluxlike[k] = np.zeros(new_shape)
                if isinstance(self.fluxlike[k], u.Quantity):
                    new.fluxlike[k] *= self.fluxlike[k].unit

            # store the binned array in the appropriate place
            if k == "uncertainty":
                # uncertainties are usually standard error on the mean
                new.fluxlike[k][:, t] = binned["uncertainty"]
            else:
                # note: all quantities are weighted the same as flux (probably inversevariance)
                new.fluxlike[k][:, t] = binned["y"]

    if (new.nwave == 0) or (new.ntime == 0):
        message = f"""
        You tried to bin {self} to {new}.

        After accounting for `minimum_acceptable_ok > {minimum_acceptable_ok}`,
        all new bins would end up with no usable data points.
        Please (a) make sure your input `Rainbow` has at least
        one wavelength and time, (b) check `.ok` accurately expresses
        which data you think are usable, (c) change the `minimum_acceptable_ok`
        keyword for `.bin` to a smaller value, and/or (d) try larger bins.
        """
        cheerfully_suggest(message)
        raise RuntimeError("No good data to bin! (see above)")

    # make sure dictionaries are on the up and up
    new._validate_core_dictionaries()

    # figure out the scales, after binning
    new._guess_wscale()
    new._guess_tscale()

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # deal with bins that are smaller than original
    N = new.wavelike["unbinned_wavelengths_per_binned_wavelength"]
    if minimum_points_per_bin is None:
        _warn_about_weird_binning(N, "wavelength")
    else:
        ok = new.wavelike.get("ok", np.ones(new.nwave, bool))
        new.wavelike["ok"] = ok * (N >= minimum_points_per_bin)

    # return the new Rainbow (with trimming if necessary)
    if trim:
        return new.trim_wavelengths(minimum_acceptable_ok=minimum_acceptable_ok)
    else:
        return new


def get_average_lightcurve_as_rainbow(self):
    """
    Produce a wavelength-integrated light curve.

    The average across wavelengths is uncertainty-weighted.

    This uses `bin`, which is a horribly slow way of doing what is
    fundamentally a very simply array calculation, because we
    don't need to deal with partial pixels.

    Returns
    -------
    lc : Rainbow
        A `Rainbow` object with just one wavelength.
    """
    h = self._create_history_entry("get_average_spectrum_as_rainbow", locals())

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        new = self.bin(nwavelengths=self.nwave, trim=False)

    new._record_history_entry(h)
    return new


def get_average_spectrum_as_rainbow(self):
    """
    Produce a time-integrated spectrum.

    The average across times is uncertainty-weighted.

    This uses `bin`, which is a horribly slow way of doing what is
    fundamentally a very simply array calculation, because we
    don't need to deal with partial pixels.

    Returns
    -------
    lc : Rainbow
        A `Rainbow` object with just one time.
    """
    h = self._create_history_entry("get_average_spectrum_as_rainbow", locals())

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        new = self.bin(ntimes=self.ntime, trim=False)

    new._record_history_entry(h)
    return new
