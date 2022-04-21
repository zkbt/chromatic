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

    w = rainbow.fluxlike["wavelength"]
    w[rainbow.ok == False] = np.nan
    dw_per_time = np.gradient(w, axis=rainbow.waveaxis)
    R_per_time = w / dw_per_time

    if visualize:
        fi, ax = plt.subplots(1, 2, figsize=(8, 3), dpi=300)
        plt.sca(ax[0])
        plt.imshow(dw_per_time, aspect="auto", vmin=0)
        plt.xlabel("Time Index")
        plt.ylabel("Wavelength Index")
        plt.title("$\Delta\lambda$")
        plt.colorbar(orientation="horizontal", pad=0.25)

        plt.sca(ax[1])
        plt.imshow(R_per_time, aspect="auto", vmin=0)
        plt.xlabel("Time Index")
        plt.ylabel("Wavelength Index")
        plt.title("R = $\lambda/\Delta\lambda$")
        plt.colorbar(orientation="horizontal", pad=0.25)
        plt.tight_layout()

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
        )
        plt.title("$\Delta\lambda$")
        plt.xlabel(f"Wavelength ({w.unit})")
        plt.ylabel("$\Delta\lambda$")
        plt.legend()

        plt.sca(ax[1])
        plt.plot(w, R_per_time, alpha=0.2)
        plt.plot(
            shared_w,
            shared_R / supersampling,
            color="black",
            label=f"(shared)/{supersampling}",
        )
        plt.title("R = $\lambda/\Delta\lambda$")
        plt.xlabel(f"Wavelength ({w.unit})")
        plt.ylabel("R = $\lambda/\Delta\lambda$")
        plt.legend()
        plt.tight_layout()
    return shared_w


def align_wavelengths(self, **kw):
    """
    Use 2D wavelength information to align onto a single 1D wavelength array.

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

    """
    # create a shared wavelength array
    shared_wavelengths = self._create_shared_wavelength_axis(**kw)

    # bin the rainbow onto that new grid, starting from 2D wavelengths
    shifted = self.bin(wavelength=shared_wavelengths, starting_wavelengths="2D")

    return shifted
