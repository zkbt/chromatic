from ..rainbows import *
from ..rainbows.actions.align_wavelengths import _create_shared_wavelength_axis
from .setup_tests import *


def create_simulation_with_wobbly_wavelengths(
    fractional_shift=0.002,
    signal_to_noise=10,
    dw=0.0001 * u.micron,
    wlim=[0.9, 1.0] * u.micron,
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
    r = SimulatedRainbow(dw=dw, wlim=wlim, **kw).inject_noise(
        signal_to_noise=signal_to_noise
    )
    wobbly_wavelengths = r.wavelength[:, np.newaxis] * np.random.normal(
        1, fractional_shift, r.ntime
    )
    r.fluxlike["wavelength"] = wobbly_wavelengths
    r.fluxlike["model"] = fake_spectrum(r.fluxlike["wavelength"])
    r.fluxlike["flux"] = r.fluxlike["flux"] * r.fluxlike["model"]

    return r


def test_create_shared_wavelength_axis(fractional_shift=0.002, dw=0.0001 * u.micron):
    r = create_simulation_with_wobbly_wavelengths(
        fractional_shift=fractional_shift, dw=dw
    )
    _create_shared_wavelength_axis(r, wscale="linear", visualize=True)
    plt.savefig(
        os.path.join(
            test_directory,
            "create-shared-wavelength-demonstration.pdf",
        )
    )


def test_align_wavelengths(fractional_shift=0.002, dw=0.0001 * u.micron):
    r = create_simulation_with_wobbly_wavelengths(
        fractional_shift=fractional_shift, dw=dw
    )
    a = r.align_wavelengths()

    fi, ax = plt.subplots(2, 2, dpi=300, figsize=(8, 6))
    for i, x in enumerate([r, a]):
        ax[0, i].imshow(x.flux)
        plt.title(["original", "aligned"][i] + " flux")
        ax[1, i].imshow(x.fluxlike["wavelength"])
        plt.title(["original", "aligned"][i] + " wavelength")
    plt.tight_layout()

    plt.savefig(os.path.join(test_directory, "wavelength-alignment-demonstration.pdf"))
    plt.close("all")


def test_align_wavelengths_with_not_ok_data(visualize=False):
    for ok_fraction in [0.25, 0.5, 0.75, 1.0]:
        r = create_simulation_with_wobbly_wavelengths(
            fractional_shift=0.0005,
            wlim=[0.99, 1.01] * u.micron,
            dw=0.001 * u.micron,
            dt=20 * u.minute,
        )
        r.fluxlike["ok"] = np.random.uniform(size=r.shape) < ok_fraction
        cautious = r.align_wavelengths(ok_threshold=1)
        carefree = r.align_wavelengths(ok_threshold=0)
        if visualize:
            cautious.imshow_quantities()
            plt.suptitle(ok_fraction)
            carefree.imshow_quantities()
            plt.suptitle(ok_fraction)

        assert np.all((cautious.ok == 1) | (cautious.ok == 0))

        if np.any(r.ok == 0):
            assert np.any((carefree.ok != 1) & (carefree.ok != 0))
