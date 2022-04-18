"""
Tools for resampling array from one grid
of independent variables to another.
"""

from .imports import *

__all__ = ["bintoR", "bintogrid", "resample_while_conserving_flux"]


def calculate_bin_widths(x):
    """
    If x is an array of bin centers, calculate the bin sizes.
    (assumes outermost bins are same size as their neighbors)

    Parameters
    ----------
    x : np.array
        The array of bin centers.

    Returns
    ----------
    s : np.array
        The array of bin sizes (total size, from left to right).
    """

    binsize = np.zeros_like(x)
    binsize[0:-1] = x[1:] - x[0:-1]
    binsize[-1] = binsize[-2]
    return binsize


def calculate_bin_edges(x):
    """
    If x is an array of bin centers, calculate the bin edges.
    (assumes outermost bins are same size as their neighbors)

    Parameters
    ----------
    x : np.array
        The array of bin centers.

    Returns
    ----------
    l : np.array
        The left edges of the bins.
    r : np.array
        The right edges of the bins.
    """

    # TODO: confirm this isn't doing something weird on log grids!

    # what are bin edges (making a guess for those on the ends)
    xbinsize = calculate_bin_widths(x)
    xleft = x - xbinsize / 2.0
    xright = x + xbinsize / 2.0

    return xleft, xright


def plot_as_boxes(x, y, **kwargs):
    """
    Plot with boxes, to show the left and right edges of a box.
    This is useful, or example, to plot flux associated with
    pixels, in case you are trying to do a sub-pixel resample
    or interpolation or shift.

    Parameters
    ----------
    x : np.array
        The original independent variable.
    y : np.array
        The original dependent variable (same size as x).
    **kwargs : dict
        All additional keywords will be passed to plt.plot
    """

    # what are bin edges (making a guess for those on the ends)
    xleft, xright = calculate_bin_edges(x)

    # create a array doubling up the y values and interleaving the edges
    plot_x = np.vstack((xleft, xright)).reshape((-1,), order="F")
    plot_y = np.vstack((y, y)).reshape((-1,), order="F")

    # plot those constructed arrays
    plt.plot(plot_x, plot_y, **kwargs)


def resample_while_conserving_flux(
    xin=None,
    yin=None,
    xout=None,
    xin_edges=None,
    xout_edges=None,
    replace_nans=0.0,
    visualize=False,
    pause=False,
):
    """
    Starting from some initial x and y, resample onto a
    different grid (either higher or lower resolution),
    while conserving total flux.

    When including the entire range of xin,
    sum(yout) == sum(yin) should be true.

    When including only part of the range of xin,
    the integral between any two points should be conserved.

    Parameters
    ----------

    xin : np.array
        The original independent variable.
    yin : np.array
        The original dependent variable (same size as x).
    xout : np.array
        The new grid of independent variables onto which
        you want to resample the y values. Refers to the
        center of each bin (use `xout_edges` for finer
        control over the exact edges of the bins)
    xin_edges : np.array
        The edges of the original independent variable bins.
        The left and right edges of the bins are interpreted
        to be `xin_edges[:-1]` and `xin_edges[1:]`,
        respectively, so the associated `yin` should have exactly
        1 fewer element than `xin_edges`. This provides finer
        control over the size of each bin in the input than
        simply supplying `xin`(still a little experimental)
        They should probably be sorted?
    xout_edges : np.array
        The edges of the new grid of bins for the independent
        variable, onto which you want to resample the y
        values. The left and right edges of the bins will be,
        respectively, `xout_edges[:-1]` and `xout_edges[1:]`,
        so the size of the output array will be
        `len(xout_edges) - 1`
    replace_nans : float, str
        Replace nan values with this value.
        `replace_nans = 0`
            will add no flux where nans are
        `replace_nans = nan`
            will ensure you get nans returned everywhere
            if you try to resample over any nan
        `replace_nans = 'interpolate'`
            will try to replace nans by linearly interpolating
            from nearby values (not yet implemented)
    visualize : bool
        Should we make a plot showing whether it worked?
    pause : bool
        Should we pause to wait for a key press?

    Returns
    -------
    result : dict
        A dictionary containing...
            `x` = the center of the output grid
            `y` = the resampled value on the output grid
            `edges` = the edges of the output grid, which will
                have one more element than x or y
    """

    # make sure there are some reasonable input options
    assert (xin is not None) or (xin_edges is not None)
    assert yin is not None
    assert (xout is not None) or (xout_edges is not None)

    # set up the bins, to calculate cumulative distribution of y
    if xin_edges is None:
        # make sure the sizes match up
        assert len(xin) == len(yin)
        # sort to make sure x is strictly increasing
        s = np.argsort(xin)
        xin_sorted = xin[s]
        yin_sorted = yin[s]
        # estimate some bin edges (might fail for non-uniform grids)
        xin_left, xin_right = calculate_bin_edges(xin_sorted)
        # define an array of edges
        xin_edges = np.hstack([xin_left, xin_right[-1]])
    else:
        # make sure the sizes match up
        assert len(xin_edges) == (len(yin) + 1)
        # sort to make sure x is strictly increasing
        s = np.argsort(xin_edges)
        xin_left, xin_right = xin_edges[s][:-1], xin_edges[s][1:]
        xin_sorted = (xin_left + xin_right) / 2
        yin_sorted = yin[s[:-1]]

    # the first element should be the left edge of the first pixel
    # last element will be right edge of last pixel
    xin_for_cdf = xin_edges

    # to the left of the first pixel, assume flux is zero
    yin_for_cdf = np.hstack([0, yin_sorted])

    # correct for any non-finite values
    bad = np.isnan(yin_for_cdf)
    if replace_nans == "interpolate":
        raise NotImplementedError(
            "The `replace_nans='interpolate'`` option doens't exist yet!"
        )
    yin_for_cdf[bad] = replace_nans

    # calculate the CDF of the flux (at pixel edge locations)
    cdfin = np.cumsum(yin_for_cdf)

    # create an interpolator for that CDF
    cdfinterpolator = interp1d(
        xin_for_cdf,
        cdfin,
        kind="linear",
        bounds_error=False,
        fill_value=(0.0, np.sum(yin)),
    )

    # calculate bin edges (of size len(xout)+1)
    if xout_edges is None:
        xout_left, xout_right = calculate_bin_edges(xout)
        xout_edges = np.hstack([xout_left, xout_right[-1]])
    else:
        xout_left, xout_right = xout_edges[:-1], xout_edges[1:]
        xout = (xout_left + xout_right) / 2

    xout_for_cdf = np.hstack([xout_left, xout_right[-1]])

    # interpolate the CDF onto those bin edges
    cdfout = cdfinterpolator(xout_for_cdf)

    # take  derivative of the CDF to get flux per resampled bin
    # (xout is bin center, and yout is the flux in that bin)
    yout = np.diff(cdfout)

    if visualize:
        fi, (ax_cdf, ax_pdf) = plt.subplots(2, 1, sharex=True, dpi=300, figsize=(8, 8))
        inkw = dict(
            color="black",
            alpha=1,
            linewidth=3,
            marker=".",
            markeredgecolor="none",
        )
        outkw = dict(
            color="darkorange",
            alpha=1,
            linewidth=1,
            marker=".",
            markersize=8,
            markeredgecolor="none",
        )

        legkw = dict(
            frameon=False,
            loc="upper left",
        )

        xinbinsize = xin_right - xin_left
        xoutbinsize = xout_right - xout_left
        # plot the PDFs
        plt.sca(ax_pdf)
        plt.ylabel("Flux per (Original) Pixel")
        plt.xlabel("Pixel")
        # plot the original pixels (in df/dpixel to compare with resampled)
        plot_as_boxes(
            xin_sorted, yin_sorted / xinbinsize, label="Original Pixels", **inkw
        )

        # what would a bad interpolation look like?
        interpolate_badly = interp1d(
            xin_sorted,
            yin_sorted / xinbinsize,
            kind="linear",
            bounds_error=False,
            fill_value=0.0,
        )
        plt.plot(
            xout,
            interpolate_badly(xout),
            color="cornflowerblue",
            alpha=1,
            linewidth=1,
            marker=".",
            markersize=8,
            markeredgecolor="none",
            label="Silly Simple Interpolation",
        )

        # plot the flux-conserving resampled data (again, in df/d"pixel")
        plt.plot(
            xout, yout / xoutbinsize, label="Flux-Conserving Interpolation", **outkw
        )

        plt.legend(**legkw)

        # plot the CDFs
        plt.sca(ax_cdf)
        plt.ylabel("Cumulative Flux (from left)")

        # plot the original CDF
        plt.plot(xin_for_cdf, cdfin, label="Original Pixels", **inkw)

        # plot the interpolated CDF
        plt.plot(xout_for_cdf, cdfout, label="Flux-Conserved Resample", **outkw)
        if pause:
            a = input(
                "Pausing a moment to check on interpolation; press return to continue."
            )

        print("{:>6} = {:.5f}".format("Actual", np.sum(yin)))
        print(
            "{:>6} = {:.5f}".format(
                "Silly",
                np.sum(interpolate_badly(xout) * xoutbinsize),
            )
        )
        print("{:>6} = {:.5f}".format("CDF", np.sum(yout)))

    # return the resampled y-values
    return {"x": xout, "edges": xout_edges, "y": yout}


def bintogrid(
    x,
    y,
    unc=None,
    newx=None,
    newx_edges=None,
    dx=None,
    weighting="inversevariance",
    drop_nans=True,
):
    """
    Bin any x and y array onto a linearly uniform grid.

    Parameters
    ----------

    x : np.array
        The original independent variable.
        (For a spectrum example = wavelength)
    y : np.array
        The original dependent variable (same size as x).
        (For a spectrum example = flux)
    unc : np.array, or None
        The unceratinty on the dependent variable
        (For a spectrum example = the flux uncertainty)
    dx : np.array
        The fixed spacing for creating a new, linearly uniform
        grid that start at the first value of x. This will
        be ignored if `newx` != None.
    newx : np.array
        A new custom grid onto which we should bin.

    weighting : str
        How should we weight values when averaging
        them together into one larger bin?
        `weighting = 'inversevariance'`
            weights = 1/unc**2
         `weighting = {literally anything else}`
            uniform weights
        This will have no impact if `unc == None`, or for any
        new bins that effectively overlap less than one original
        unbinned point.
    drop_nans : bool
        Should we skip any bins turn out to be nans?
        This most often happens when bins are empty.

    Returns
    -------
    result : dict
        A dictionary containing at least...
            `x` = the center of the output grid
            `y` = the resampled value on the output grid
            `edges` = the edges of the output grid, which will
                have one more element than x or y
        ...and possibly also
            `uncertainty` = the calculated uncertainty per bin

    """

    try:
        x_unit = x.unit
        x_without_unit = x.value
    except AttributeError:
        x_unit = 1
        x_without_unit = x

    try:
        y_unit = y.unit
        y_without_unit = y.value
    except AttributeError:
        y_unit = 1
        y_without_unit = y

    # make up a grid, if one wasn't specified
    if newx is None:
        dx_without_unit = u.Quantity(dx).to(x_unit).value
        newx_without_unit = np.arange(
            np.nanmin(x_without_unit),
            np.nanmax(x_without_unit) + dx_without_unit,
            dx_without_unit,
        )
    else:
        newx_without_unit = u.Quantity(newx).to(x_unit).value

    # don't complain about zero-divisions in here (to allow infinite uncertainties)
    with np.errstate(divide="ignore", invalid="ignore"):

        # calculate weight integrals for the bin array
        ok = np.isnan(y_without_unit) == False

        # resample the sums onto that new grid
        if unc is None:
            weights = np.ones_like(x_without_unit)
        else:
            if weighting == "inversevariance":
                weights = 1 / unc ** 2
            else:
                weights = np.ones_like(x_without_unit)

            # ignore infinite weights (= 0 uncertainties)
            ok *= np.isfinite(weights)

        if np.any(ok):
            # TO-DO: check this nan handling on input arrays is OK?
            numerator = resample_while_conserving_flux(
                x_without_unit[ok], (y_without_unit * weights)[ok], newx_without_unit
            )
            denominator = resample_while_conserving_flux(
                x_without_unit[ok], (weights)[ok], newx_without_unit
            )

            # the binned weighted means on the new grid
            newy = numerator["y"] / denominator["y"]

            # the standard error on the means, for those bins
            newunc = np.sqrt(1 / denominator["y"])
        else:
            newy = np.nan * newx_without_unit
            newunc = np.nan * newx_without_unit

    # remove any empty bins
    if drop_nans:
        ok = np.isfinite(newy)
    else:
        ok = np.ones_like(newx_without_unit).astype(bool)

    # if no uncertainties were given, don't return uncertainties
    result = {}
    result["x"] = newx_without_unit[ok] * x_unit
    result["y"] = newy[ok] * y_unit
    result["edges"] = numerator["edges"] * x_unit
    if unc is not None:
        result["uncertainty"] = newunc[ok] * y_unit
    return result


def bintoR(
    x, y, unc=None, R=50, xlim=None, weighting="inversevariance", drop_nans=True
):
    """
    Bin any x and y array onto a logarithmicly uniform grid,
    characterized by

    Parameters
    ----------

    x : np.array
        The original independent variable.
        (For a spectrum example = wavelength)
    y : np.array
        The original dependent variable (same size as x).
        (For a spectrum example = flux)
    unc : np.array, or None
        The unceratinty on the dependent variable
        (For a spectrum example = the flux uncertainty)
    R : np.array
        The spectral resolution R=x/dx for creating a new,
        logarithmically uniform grid that starts at the first
        value of x.
    xlim : list, np.array
        A two-element list indicating the min and max values of
        x for the new logarithmically spaced grid. If None,
        these limits will be created from the data themselves
    weighting : str
        How should we weight values when averaging
        them together into one larger bin?
        `weighting = 'inversevariance'`
            weights = 1/unc**2
         `weighting = {literally anything else}`
            uniform weights
        This will have no impact if `unc == None`, or for any
        new bins that effectively overlap less than one original
        unbinned point.
    drop_nans : bool
        Should we skip any bins turn out to be nans?
        This most often happens when bins are empty.

    Returns
    -------
    newx : np.array
        The new grid of x values.
    newy : np.array
        The new gridded y values.
    newunc : np.array
        The new gridded y uncertainties.
        These will be returned only if `unc != None` in the
        original inputs. Otherwise, only two outputs will be
        returned, newx and newy.

    # TODO: confirm nans in input arrays are handled OK
    """

    try:
        x_unit = x.unit
        x_without_unit = x.value
    except AttributeError:
        x_unit = 1
        x_without_unit = x

    # create a new grid of x at the given resolution
    lnx = np.log(x_without_unit)
    dnewlnx = 1.0 / R

    # set the limits of the new xgrid (in log space)
    if xlim is None:
        # use the input grid to set the limits
        lnxbottom, lnxtop = np.nanmin(lnx), np.nanmax(lnx)
    else:
        # use the custom xlim to set the limits
        lnxbottom, lnxtop = xlim

    # create a new, log-uniform grid of x values
    newlnx = np.arange(lnxbottom, lnxtop + dnewlnx, dnewlnx)

    # now do the binning on a uniform grid of lnx
    result = bintogrid(
        lnx, y, unc, newx=newlnx, weighting=weighting, drop_nans=drop_nans
    )

    # convert back from log to real values
    result["x"] = np.exp(result["x"]) * x_unit
    result["edges"] = np.exp(result["edges"]) * x_unit

    return result
