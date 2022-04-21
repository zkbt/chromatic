from ..rainbows import *
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
            f *= 1 - a * np.exp(-0.5 * (w - c) ** 2 / s ** 2)
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
