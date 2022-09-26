from ..rainbows import *
from ..rainbows.actions.align_wavelengths import _create_shared_wavelength_axis
from .setup_tests import *


def create_simulation_with_wobbly_wavelengths(
    fractional_shift=0.002,
    signal_to_noise=100,
    dw=0.001 * u.micron,
    wlim=[0.95, 1.05] * u.micron,
    dt=30 * u.minute,
    **kw
):

    # set up a function to create a fake absorption line spectrum
    N = 10
    centers = (
        np.random.uniform(wlim[0].to_value("micron"), wlim[1].to_value("micron"), N)
        * u.micron
    )
    sigmas = np.random.uniform(0.001, 0.003, N) * u.micron
    amplitudes = np.random.uniform(0, 0.5, N)

    def fake_spectrum(w):
        f = u.Quantity(np.ones(np.shape(w)))
        for c, s, a in zip(centers, sigmas, amplitudes):
            f *= 1 - a * np.exp(-0.5 * (w - c) ** 2 / s**2)
        return f

    # create
    r = SimulatedRainbow(dw=dw, wlim=wlim, dt=dt, **kw).inject_noise(
        signal_to_noise=signal_to_noise
    )
    wobbly_wavelengths = r.wavelength[:, np.newaxis] * np.random.normal(
        1, fractional_shift, r.ntime
    )
    r.fluxlike["wavelength_2d"] = wobbly_wavelengths
    r.fluxlike["model"] = fake_spectrum(r.fluxlike["wavelength_2d"])
    r.fluxlike["flux"] = r.fluxlike["flux"] * r.fluxlike["model"]

    return r


def test_create_shared_wavelength_axis(fractional_shift=0.002, dw=0.01 * u.micron):
    r = create_simulation_with_wobbly_wavelengths(
        fractional_shift=fractional_shift, dw=dw
    )
    _create_shared_wavelength_axis(r, wscale="linear", visualize=True)
    plt.savefig(
        os.path.join(
            test_directory,
            "demonstration-of-creating-shared-wavelength-axis.pdf",
        )
    )


def test_align_wavelengths(fractional_shift=0.002, dw=0.001 * u.micron, **kw):
    r = create_simulation_with_wobbly_wavelengths(
        fractional_shift=fractional_shift, dw=dw, **kw
    )
    a = r.align_wavelengths()
    plt.close("all")
    fi, ax = plt.subplots(2, 2, dpi=300, figsize=(8, 6), constrained_layout=True)
    for i, x in enumerate([r, a]):
        plt.sca(ax[0, i])
        plt.imshow(x.flux, aspect="auto", cmap="gray")
        plt.title(["original", "aligned"][i] + " flux")
        plt.xlabel("time index")
        plt.ylabel("wavelength index")
        plt.colorbar()

        plt.sca(ax[1, i])
        plt.imshow(x.fluxlike["wavelength_2d"], aspect="auto")
        plt.title(["original", "aligned"][i] + " wavelength")
        plt.xlabel("time index")
        plt.ylabel("wavelength index")
        plt.colorbar()

    plt.savefig(
        os.path.join(
            test_directory,
            "demonstration-of-aligning-2D-wavelengths-to-shared-axis.pdf",
        )
    )


def test_align_wavelengths_with_not_ok_data(visualize=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for ok_fraction in [0.25, 0.5, 0.75, 1.0]:
            r = create_simulation_with_wobbly_wavelengths(
                fractional_shift=0.0005,
                wlim=[0.99, 1.01] * u.micron,
                dw=0.001 * u.micron,
                dt=20 * u.minute,
            )
            r.fluxlike["ok"] = np.random.uniform(size=r.shape) < ok_fraction
            cautious = r.align_wavelengths(minimum_acceptable_ok=1)
            carefree = r.align_wavelengths(minimum_acceptable_ok=0)
            if visualize:
                cautious.imshow_quantities()
                plt.suptitle(ok_fraction)
                carefree.imshow_quantities()
                plt.suptitle(ok_fraction)

        # if np.any(r.ok == 0):
        #    assert np.any((carefree.ok != 1) & (carefree.ok != 0))
