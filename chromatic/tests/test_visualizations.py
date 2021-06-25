from ..rainbows import *
from .setup_tests import *


def test_imshow():
    plt.figure()
    SimulatedRainbow(R=10).imshow()

    plt.figure()
    SimulatedRainbow(dw=0.2 * u.micron).imshow()

    plt.figure()
    Rainbow(wavelength=np.arange(1,5)*u.micron,
            time=np.arange(1,6)*u.hour,
            flux=np.ones((4,5))).imshow()

    fi, ax = plt.subplots(2, 1, sharex=True)
    SimulatedRainbow(R=10).imshow(w_unit="nm", ax=ax[0])
    SimulatedRainbow(dw=0.2 * u.micron).imshow(ax=ax[1], w_unit="nm")
    plt.savefig(os.path.join(test_directory, "imshow-demonstration.pdf"))


def test_plot():
    SimulatedRainbow(R=10).plot()
    plt.savefig(os.path.join(test_directory, "plot-demonstration.pdf"))
