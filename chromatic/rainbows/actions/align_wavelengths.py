from ...imports import *
from ...resampling import *

__all__ = ["_create_shared_wavelength_axis", "align_wavelengths"]


def _create_shared_wavelength_axis(
    rainbow, wscale="linear", supersampling=1, visualize=False
):
    """
    Create a shared 1D wavelength axis that tries to compress
    different wavelength axes associated with different times.

    Parameters
    ----------
    wscale : str
        What kind of a new wavelength axis should be created?
        Options include:
            'linear' = constant d[wavelength] between grid points
            'log' = constant d[wavelength]/[wavelength] between grid points

    supersampling : float
        By how many times should we increase or decrease the wavelength sampling?
        In general, values >1 will split each input wavelength grid point into
        multiple supersampled wavelength grid points, values close to 1 will
        produce approximately one output wavelength for each input wavelength,
        and values <1 will average multiple input wavelengths into a single output
        wavelength bin.

        Unless this is significantly less than 1, there's a good chance your output
        array may have strong correlations between one or more adjacent wavelengths.
        Be careful when trying to use the resulting uncertainties!

        (FIXME = add a way to estimate a covariance matrix when binning?)

    visualize : bool
        Should we make some plots showing how the shared wavelength
        axis compares to the original input wavelength axes?
    """
    w = rainbow.fluxlike["wavelength_2d"] * 1
    w[rainbow.ok == False] = np.nan
    dw_per_time = np.gradient(w, axis=rainbow.waveaxis)
    R_per_time = w / dw_per_time

    if visualize:
        fi, ax = plt.subplots(1, 2, figsize=(8, 3), dpi=300)
        plt.sca(ax[0])
        plt.imshow(dw_per_time, aspect="auto", vmin=0)
        plt.xlabel("Time Index")
        plt.ylabel("Wavelength Index")
        plt.title(r"$\Delta\lambda$")
        plt.colorbar(orientation="horizontal", pad=0.25)

        plt.sca(ax[1])
        plt.imshow(R_per_time, aspect="auto", vmin=0)
        plt.xlabel("Time Index")
        plt.ylabel("Wavelength Index")
        plt.title(r"R = $\lambda/\Delta\lambda$")
        plt.colorbar(orientation="horizontal", pad=0.25)

    min_w, max_w = np.nanmin(w).to("micron").value, np.nanmax(w).to("micron").value
    if wscale == "linear":
        dw = np.nanmedian(dw_per_time).to("micron").value / supersampling
        shared_w = np.arange(min_w, max_w + dw, dw) * u.micron
    elif wscale == "log":
        R = np.nanmedian(w / dw_per_time) * supersampling
        shared_w = np.exp(np.arange(np.log(min_w), np.log(max_w) + 1 / R, 1 / R))
    shared_dw = np.gradient(shared_w)
    shared_R = shared_w / shared_dw

    if visualize:
        fi, ax = plt.subplots(1, 2, figsize=(8, 3), dpi=300)
        plt.sca(ax[0])
        plt.plot(w, dw_per_time, alpha=0.2)
        plt.plot(
            shared_w,
            shared_dw * supersampling,
            color="black",
            label=f"{supersampling}x(shared)",
            marker=".",
        )
        plt.title(r"$\Delta\lambda$")
        plt.xlabel(f"Wavelength ({w.unit})")
        plt.ylabel(r"$\Delta\lambda$")
        plt.legend()

        plt.sca(ax[1])
        plt.plot(w, R_per_time, alpha=0.2)
        plt.plot(
            shared_w,
            shared_R / supersampling,
            color="black",
            label=f"(shared)/{supersampling}",
            marker=".",
        )
        plt.title(r"R = $\lambda/\Delta\lambda$")
        plt.xlabel(f"Wavelength ({w.unit})")
        plt.ylabel(r"R = $\lambda/\Delta\lambda$")
        plt.legend()
    return shared_w


def align_wavelengths(
    self,
    minimum_acceptable_ok=1,
    minimum_points_per_bin=0,
    wscale="linear",
    supersampling=1,
    visualize=False,
):
    """
    Use 2D wavelength information to align onto a single 1D wavelength array.

    This relies on the existence of a `.fluxlike['wavelength_2d']` array,
    expressing the wavelength associated with each flux element.
    Those wavelengths will be used to (a) establish a new compromise
    wavelength grid and (b) bin the individual timepoints onto that
    new grid, effectively shifting the wavelengths to align.

    Parameters
    ----------
    minimum_acceptable_ok : float, optional
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
                  will go into the binning. All other points
                  will effectively be interpolated over. Flux
                  uncertainties *should* be inflated appropriately,
                  but it's very possible to create correlated
                  bins next to each other if many of your ingoing
                  data points are not perfectly OK.
            minimum_acceptable_ok = 1
                  All data points that aren't definitely not OK
                  will go into the binning. The OK-ness of points
                  will propagate onward for future binning.
            minimum_acceptable_ok < 0
                  All data points will be included in the bin.
                  The OK-ness will propagate onward.
    wscale : str, optional
        What kind of a new wavelength axis should be created?
        Options include:
            'linear' = constant d[wavelength] between grid points
            'log' = constant d[wavelength]/[wavelength] between grid points
    supersampling : float, optional
        By how many times should we increase or decrease the wavelength sampling?
        In general, values >1 will split each input wavelength grid point into
        multiple supersampled wavelength grid points, values close to 1 will
        produce approximately one output wavelength for each input wavelength,
        and values <1 will average multiple input wavelengths into a single output
        wavelength bin.
        Unless this is significantly less than 1, there's a good chance your output
        array may have strong correlations between one or more adjacent wavelengths.
        Be careful when trying to use the resulting uncertainties!
    visualize : bool
        Should we make some plots showing how the shared wavelength
        axis compares to the original input wavelength axes?

    Returns
    -------
    rainbow : RainbowWithModel
        A new `RainbowWithModel` object, with the model attached.
    """
    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("align_wavelengths", locals())

    if "wavelength_2d" not in self.fluxlike:
        cheerfully_suggest(
            f"""
        No 2D wavelength information was found, so
        it's assumed wavelengths don't need to be aligned.
        Wavelength alignment is being skipped!
        """
        )
        shifted = self._create_copy()
    else:
        # create a shared wavelength array
        shared_wavelengths = self._create_shared_wavelength_axis(
            wscale=wscale, supersampling=supersampling, visualize=visualize
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # bin the rainbow onto that new grid, starting from 2D wavelengths
            shifted = self.bin_in_wavelength(
                wavelength=shared_wavelengths,
                minimum_acceptable_ok=minimum_acceptable_ok,
                starting_wavelengths="2D",
                minimum_points_per_bin=minimum_points_per_bin,
            )

    # append the history entry to the new Rainbow
    shifted._record_history_entry(h)

    # return the new Rainbow
    return shifted
