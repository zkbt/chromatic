"""
Tools for resampling array from one grid
of independent variables to another.
"""

from .imports import *

__all__ = [
    "bintoR",
    "bintogrid",
    "resample_while_conserving_flux",
    "leftright_to_edges",
    "edges_to_leftright",
    "calculate_bin_leftright",
]


def calculate_bin_leftright(x):
    """
    If x is an array of bin centers, calculate the bin edges.
    (assumes outermost bins are same size as their neighbors)

    Parameters
    ----------
    x : array
        The array of bin centers.

    Returns
    ----------
    l : array
        The left edges of the bins.
    r : array
        The right edges of the bins.
    """

    # what are bin edges (making a guess for those on the ends)
    # xbinsize = calculate_bin_widths(x)
    # left = x - xbinsize / 2.0
    # right = x + xbinsize / 2.0

    # weird corner case!
    if len(x) == 1:
        left, right = np.sort([0, 2 * x[0]])
        return np.array([left]), np.array([right])

    inner_edges = 0.5 * np.diff(x) + x[:-1]
    first_edge = x[0] - (inner_edges[0] - x[0])
    last_edge = x[-1] + (x[-1] - inner_edges[-1])

    left = np.hstack([first_edge, inner_edges])
    right = np.hstack([inner_edges, last_edge])

    return left, right


def calculate_bin_widths(x):
    """
    If x is an array of bin centers, calculate the bin sizes.
    (assumes outermost bins are same size as their neighbors)

    Parameters
    ----------
    x : array
        The array of bin centers.

    Returns
    ----------
    s : array
        The array of bin sizes (total size, from left to right).
    """

    # OLD VERSION
    # binsize = np.zeros_like(x)
    # binsize[0:-1] = x[1:] - x[0:-1]
    # binsize[-1] = binsize[-2]
    left, right = calculate_bin_leftright(x)
    binsize = right - left
    return binsize


def plot_as_boxes(x, y, xleft=None, xright=None, **kwargs):
    """
    Plot with boxes, to show the left and right edges of a box.
    This is useful, or example, to plot flux associated with
    pixels, in case you are trying to do a sub-pixel resample
    or interpolation or shift.

    Parameters
    ----------
    x : array
        The original independent variable.
    y : array
        The original dependent variable (same size as x).
    **kwargs : dict
        All additional keywords will be passed to plt.plot
    """

    # what are bin edges (making a guess for those on the ends)
    if (xleft is None) and (xright is None):
        xleft, xright = calculate_bin_leftright(x)

    # create a array doubling up the y values and interleaving the edges
    plot_x = np.vstack((xleft, xright)).reshape((-1,), order="F")
    plot_y = np.vstack((y, y)).reshape((-1,), order="F")

    # plot those constructed arrays
    plt.plot(plot_x, plot_y, **kwargs)


def edges_to_leftright(edges):
    """
    Convert N+1 contiguous edges to two arrays of N left/right edges.
    """
    left, right = edges[:-1], edges[1:]
    return left, right


def leftright_to_edges(left, right):
    """
    Convert two arrays of N left/right edges to N+1 continugous edges.
    """
    edges = np.hstack([left, right[-1]])
    return edges


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

    When including the entire range of `xin`,
    `sum(yout) == sum(yin)` should be true.

    When including only part of the range of `xin`,
    the integral between any two points should be conserved.

    Parameters
    ----------
    xin : array
        The original independent variable.
    yin : array
        The original dependent variable (same size as x).
    xout : array
        The new grid of independent variables onto which
        you want to resample the y values. Refers to the
        center of each bin (use `xout_edges` for finer
        control over the exact edges of the bins)
    xin_edges : array
        The edges of the original independent variable bins.
        The left and right edges of the bins are interpreted
        to be `xin_edges[:-1]` and `xin_edges[1:]`,
        respectively, so the associated `yin` should have exactly
        1 fewer element than `xin_edges`. This provides finer
        control over the size of each bin in the input than
        simply supplying `xin`(still a little experimental)
        They should probably be sorted?
    xout_edges : array
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
        xin_left, xin_right = calculate_bin_leftright(xin_sorted)
        # define an array of edges
        xin_edges = leftright_to_edges(xin_left, xin_right)
    else:
        # make sure the sizes match up
        assert len(xin_edges) == (len(yin) + 1)
        # sort to make sure x is strictly increasing
        s = np.argsort(xin_edges)
        xin_left, xin_right = edges_to_leftright(xin_edges[s])
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
        xout_left, xout_right = calculate_bin_leftright(xout)
        xout_edges = leftright_to_edges(xout_left, xout_right)
    else:
        xout_left, xout_right = edges_to_leftright(xout_edges)
        xout = (xout_left + xout_right) / 2

    xout_for_cdf = leftright_to_edges(xout_left, xout_right)

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
    return {"x": xout, "x_edge_lower": xout_left, "x_edge_upper": xout_right, "y": yout}


def bintogrid(
    x=None,
    y=None,
    unc=None,
    newx=None,
    newx_edges=None,
    dx=None,
    nx=None,
    weighting="inversevariance",
    drop_nans=True,
    x_edges=None,
    visualize=False,
):
    """
    Bin any x and y array onto a linearly uniform grid.

    Parameters
    ----------
    x : array
        The original independent variable.
        (For a spectrum example = wavelength)
    y : array
        The original dependent variable (same size as x).
        (For a spectrum example = flux)
    unc : array, None
        The unceratinty on the dependent variable
        (For a spectrum example = the flux uncertainty)
    nx : array
        The number of bins from the original grid to
        bin together into the new one.
    dx : array
        The fixed spacing for creating a new, linearly uniform
        grid that start at the first value of x. This will
        be ignored if `newx` != None.
    newx : array
        A new custom grid onto which we should bin.
    newx_edges : array
        The edges of the new grid of bins for the independent
        variable, onto which you want to resample the y
        values. The left and right edges of the bins will be,
        respectively, `newx_edges[:-1]` and `newx_edges[1:]`,
        so the size of the output array will be
        `len(newx_edges) - 1`
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
    x_edges : array
        The edges of the original independent variable bins.
        The left and right edges of the bins are interpreted
        to be `x_edges[:-1]` and `x_edges[1:]`,
        respectively, so the associated `y` should have exactly
        1 fewer element than `x_edges`. This provides finer
        control over the size of each bin in the input than
        simply supplying `x`(still a little experimental)

    Returns
    -------
    result : dict
        A dictionary containing at least...
            `x` = the center of the output grid
            `y` = the resampled value on the output grid
            `x_edge_lower` = the lower edges of the output grid
            `x_edge_upper` = the upper edges of the output grid
        ...and possibly also
            `uncertainty` = the calculated uncertainty per bin


    The order of precendence for setting the new grid is
    [`newx_edges`, `newx`, `dx`, `nx`]
    The first will be used, and others will be ignored.
    """

    # check that an OK set of inputs has been supplied
    if (x is not None) and (x_edges is not None):
        raise RuntimeError(
            """🌈 Both `x` and `x_edges` were supplied to `bintogrid`. Confusing!"""
        )
    if (x is None) and (x_edges is None):
        raise RuntimeError(
            """🌈 At least one of `x` or `x_edges` must be supplied to `bintogrid`."""
        )
    if y is None:
        raise RuntimeError("""🌈 `y` must be supplied to `bintogrid`.""")

    # make sure the edges and the centers are set
    if x is None:
        x_left, x_right = edges_to_leftright(x_edges)
        x = 0.5 * (left + right)
    else:
        x_left, x_right = calculate_bin_leftright(x)
        x_edges = leftright_to_edges(x_left, x_right)
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

    # warn if multiple inputs are provided
    number_of_grid_options = np.sum([z is not None for z in [newx_edges, newx, dx, nx]])
    if number_of_grid_options > 1:
        cheerfully_suggest(
            """More than one output grid sent to `bintogrid`.
                         The one being used is the first to appear in
                         [`newx_edges`, `newx`, `dx`, `nx`]
                         but you might want to choose more carefully."""
        )

    # define inputs based on the following order
    if newx_edges is not None:
        # define grid by its edges (and define others from there)
        newx_edges_without_unit = u.Quantity(newx_edges).to(x_unit).value
        dx_without_unit = np.diff(newx_edges_without_unit)
        newx_without_unit = newx_edges_without_unit[:-1] + 0.5 * dx_without_unit
        newx_left_without_unit = newx_edges_without_unit[:-1]
        newx_right_without_unit = newx_edges_without_unit[1:]

        # make sure the final output grid is defined
        final_newx, final_newx_left, final_newx_right = (
            newx_without_unit * x_unit,
            newx_left_without_unit * x_unit,
            newx_right_without_unit * x_unit,
        )
    elif newx is not None:
        # define grid by its centers (and define others from there)
        newx_without_unit = u.Quantity(newx).to(x_unit).value
        newx_left_without_unit, newx_right_without_unit = calculate_bin_leftright(
            newx_without_unit
        )
        newx_edges_without_unit = np.hstack(
            [newx_left_without_unit, newx_right_without_unit[-1]]
        )
        dx_without_unit = np.diff(newx_edges_without_unit)

        # make sure the final output grid is defined
        final_newx, final_newx_left, final_newx_right = (
            newx_without_unit * x_unit,
            newx_left_without_unit * x_unit,
            newx_right_without_unit * x_unit,
        )
    elif dx is not None:
        # define grid by a bin width (and define others from there)
        dx_without_unit = u.Quantity(dx).to(x_unit).value
        newx_without_unit = np.arange(
            np.nanmin(x_without_unit),
            np.nanmax(x_without_unit) + dx_without_unit,
            dx_without_unit,
        )
        newx_left_without_unit, newx_right_without_unit = calculate_bin_leftright(
            newx_without_unit
        )
        newx_edges_without_unit = np.hstack(
            [newx_left_without_unit, newx_right_without_unit[-1]]
        )

        # make sure the final output grid is defined
        final_newx, final_newx_left, final_newx_right = (
            newx_without_unit * x_unit,
            newx_left_without_unit * x_unit,
            newx_right_without_unit * x_unit,
        )

    elif nx is not None:
        # keep track of the original input x values
        original_x_without_unit = x_without_unit

        # redefine the input x to indices, to do interpolation in index space
        x_without_unit = np.arange(0, len(x_without_unit))

        # define a grid of edges that will enclose the right number of indices
        x_left_i, x_right_i = calculate_bin_leftright(x_without_unit)
        newx_edges_without_unit = leftright_to_edges(x_left_i, x_right_i)[::nx]
        newx_without_unit = 0.5 * (
            newx_edges_without_unit[1:] + newx_edges_without_unit[:-1]
        )

        # calculate the actual x values corresponding to the bins
        original_edges = leftright_to_edges(
            *calculate_bin_leftright(original_x_without_unit)
        )
        final_edges = original_edges[::nx] * x_unit
        final_newx_left, final_newx_right = edges_to_leftright(final_edges)
        final_newx = 0.5 * (final_newx_left + final_newx_right)
        dx_without_unit = (final_newx_right - final_newx_left) / x_unit
    else:
        raise RuntimeError(
            """No output grid sent to `bintogrid`.
                              Please choose one of the following:
                              [`newx_edges`, `newx`, `dx`, `nx`]"""
        )

    # don't complain about zero-divisions in here (to allow infinite uncertainties)
    with np.errstate(divide="ignore", invalid="ignore"):

        # calculate weight integrals for the bin array
        ok = np.isnan(y_without_unit) == False

        # resample the sums onto that new grid
        if unc is None:
            weights = np.ones_like(x_without_unit)
        else:
            if weighting == "inversevariance":
                weights = 1 / unc**2
            else:
                weights = np.ones_like(x_without_unit)

            # ignore infinite weights (= 0 uncertainties)
            ok *= np.isfinite(weights)

        if np.any(ok):
            numerator = resample_while_conserving_flux(
                xin=x_without_unit[ok],
                yin=(y_without_unit * weights)[ok],
                xout_edges=newx_edges_without_unit,
            )
            denominator = resample_while_conserving_flux(
                xin=x_without_unit[ok],
                yin=weights[ok],
                xout_edges=newx_edges_without_unit,
            )

            # the binned weighted means on the new grid
            newy = numerator["y"] / denominator["y"]

            # the standard error on the means, for those bins
            newunc = np.sqrt(1 / denominator["y"])

            # keep track of the number of original bins going into each new bin
            number_of_original_bins_per_new_bin = resample_while_conserving_flux(
                xin=x_without_unit[ok],
                yin=np.ones_like(y_without_unit)[ok],
                xout_edges=newx_edges_without_unit,
            )["y"]
        else:
            newy = np.nan * newx_without_unit
            newunc = np.nan * newx_without_unit
            number_of_original_bins_per_new_bin = np.zeros_like(newx_without_unit)

    # remove any empty bins
    if drop_nans:
        ok = np.isfinite(newy)
    else:
        ok = np.ones_like(newx_without_unit).astype(bool)

    # if no uncertainties were given, don't return uncertainties
    result = {}

    # populate the new grid centers + edges + values
    result["x"] = final_newx[ok]
    result["x_edge_lower"] = final_newx_left[ok]
    result["x_edge_upper"] = final_newx_right[ok]

    # populate the new grid values
    result["y"] = newy[ok] * y_unit

    # populate the new grid value uncertainties
    if unc is not None:
        result["uncertainty"] = newunc[ok] * y_unit

    # store how many of the original pixels made it into this new one
    result["N_unbinned/N_binned"] = number_of_original_bins_per_new_bin[ok]
    if visualize:
        fi, ax = plt.subplots(
            2, 1, figsize=(8, 4), dpi=300, gridspec_kw=dict(height_ratios=[1, 0.2])
        )
        plt.sca(ax[0])
        plot_as_boxes(x, y, xleft=x_left, xright=x_right, color="silver", linewidth=1)
        ekw = dict(elinewidth=1, linewidth=0)
        plt.errorbar(x, y, yerr=unc, color="silver", marker="s", **ekw)
        plt.errorbar(
            result["x"],
            result["y"],
            yerr=result.get("uncertainty", None),
            xerr=0.5 * (result["x_edge_upper"] - result["x_edge_lower"]) * x_unit,
            marker="o",
            color="black",
            zorder=100,
            **ekw,
        )
        plt.sca(ax[1])
        plot_as_boxes(
            result["x"],
            result["N_unbinned/N_binned"],
            xleft=result["x_edge_lower"],
            xright=result["x_edge_upper"],
        )
        plt.ylabel("$N_{unbinned}/N_{binned}$")
        plt.ylim(0, None)

    return result


def bintoR(
    x, y, unc=None, R=50, xlim=None, weighting="inversevariance", drop_nans=True
):
    """
    Bin any x and y array onto a logarithmicly uniform grid.

    Parameters
    ----------
    x : array
        The original independent variable.
        (For a spectrum example = wavelength)
    y : array
        The original dependent variable (same size as x).
        (For a spectrum example = flux)
    unc : array, None, optional
        The unceratinty on the dependent variable
        (For a spectrum example = the flux uncertainty)
    R : array, optional
        The spectral resolution R=x/dx for creating a new,
        logarithmically uniform grid that starts at the first
        value of x.
    xlim : list, array, optional
        A two-element list indicating the min and max values of
        x for the new logarithmically spaced grid. If None,
        these limits will be created from the data themselves
    weighting : str, optional
        How should we weight values when averaging
        them together into one larger bin?
        `weighting = 'inversevariance'`
            weights = 1/unc**2
         `weighting = {literally anything else}`
            uniform weights
        This will have no impact if `unc == None`, or for any
        new bins that effectively overlap less than one original
        unbinned point.
    drop_nans : bool, optional
        Should we skip any bins turn out to be nans?
        This most often happens when bins are empty.

    Returns
    -------
    result : dict
        A dictionary containing at least...
            `x` = the center of the output grid
            `y` = the resampled value on the output grid
            `x_edge_lower` = the lower edges of the output grid
            `x_edge_upper` = the upper edges of the output grid
        ...and possibly also
            `uncertainty` = the calculated uncertainty per bin
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
    for k in ["x", "x_edge_lower", "x_edge_upper"]:
        result[k] = np.exp(result[k]) * x_unit

    return result
