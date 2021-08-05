from ..rainbows import *
from .setup_tests import *


def test_multi():
    a = SimulatedRainbow(signal_to_noise=400).inject_transit()
    b = SimulatedRainbow(signal_to_noise=800).inject_transit()

    m = MultiRainbow([a, b])
    m.bin(R=5).plot(spacing=0.02)
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
    m.align_wavelengths().wavelength
    m[:, :]
