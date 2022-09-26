from ..rainbows import *
from .setup_tests import *


def test_multi():
    a = (
        SimulatedRainbow(R=3, dt=20 * u.minute)
        .inject_transit()
        .inject_noise(signal_to_noise=400)
    )
    b = (
        SimulatedRainbow(R=3, dt=20 * u.minute)
        .inject_transit()
        .inject_noise(signal_to_noise=800)
    )

    m = MultiRainbow([a, b])
    m.plot(spacing=0.02)
    plt.savefig(os.path.join(test_directory, "multi-plot-demonstration.png"))

    m.imshow()
    plt.savefig(os.path.join(test_directory, "multi-imshow-demonstration.png"))

    m.animate_lightcurves(
        filename=os.path.join(
            test_directory, "multi-animate-lightcurves-demonstration.gif"
        )
    )
    m.animate_spectra(
        filename=os.path.join(test_directory, "multi-animate-spectra-demonstration.gif")
    )

    m.normalize()
    # m.align_wavelengths().wavelength
    # m[:, :]
    plt.close("all")


def test_compare_wrappers():
    rainbows = [SimulatedRainbow(R=10 ** np.random.uniform(0.1, 2)) for _ in range(3)]

    a = compare_rainbows(rainbows)
    b = rainbows[0].compare(rainbows)

    assert a.names == b.names
    assert a.rainbows == b.rainbows
