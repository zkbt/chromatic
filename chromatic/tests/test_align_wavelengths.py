from ..rainbows import *
from ..rainbows.actions.align_wavelengths import _create_shared_wavelength_axis
from .setup_tests import *


def create_simulation_with_wobbly_wavelengths(
    fractional_shift=0.002,
    signal_to_noise=10,
    dw=0.0001 * u.micron,
    wlim=[0.9, 1.0] * u.micron,
):

    # set up a function to create a fake absorption line spectrum
    N = 10
    centers = np.random.uniform(0.9, 1, N) * u.micron
    sigmas = np.random.uniform(0.001, 0.003, N) * u.micron
    amplitudes = np.random.uniform(0, 0.5, N)

    def fake_spectrum(w):
        f = u.Quantity(np.ones(np.shape(w)))
        for c, s, a in zip(centers, sigmas, amplitudes):
            f *= 1 - a * np.exp(-0.5 * (w - c) ** 2 / s**2)
        return f

    # create
    r = SimulatedRainbow(signal_to_noise=signal_to_noise, dw=dw, wlim=wlim)
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
