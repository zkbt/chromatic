from ..rainbows import *
from .setup_tests import *


def test_normalize(plot=False):
    N = 37
    w = np.logspace(0, 1, N) * u.micron
    f = np.cos(2 * np.pi * w.value / 3) + 1

    snr = 100
    a = SimulatedRainbow(wavelength=w).inject_noise(signal_to_noise=snr)
    b = SimulatedRainbow(wavelength=w, star_flux=f).inject_noise(signal_to_noise=snr)
    c = (
        SimulatedRainbow(wavelength=w, star_flux=f)
        .inject_transit()
        .inject_noise(signal_to_noise=snr)
    )

    for x in [a, b, c]:
        nw = x.normalize(axis="w")
        nt = x.normalize(axis="t")
        nwt = x.normalize(axis="w").normalize(axis="t")
        ntw = x.normalize(axis="t").normalize(axis="w")

        for r in [nw, nt, nwt, ntw]:
            r.fluxlike["relative-uncertainty"] = r.uncertainty / r.flux
            assert np.all(np.isclose(r.uncertainty / r.flux, 1 / snr, rtol=0.1))
            if plot:
                r.imshow_quantities(maxcol=4)
    plt.close("all")
