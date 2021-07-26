from ..rainbows import *
from .setup_tests import *


def test_imshow():
    plt.figure()
    SimulatedRainbow(R=10).imshow()

    plt.figure()
    SimulatedRainbow(dw=0.2 * u.micron).imshow()

    plt.figure()
    Rainbow(
        wavelength=np.arange(1, 5) * u.micron,
        time=np.arange(1, 6) * u.hour,
        flux=np.ones((4, 5)),
    ).imshow()

    fi, ax = plt.subplots(2, 1, sharex=True)
    SimulatedRainbow(R=10).imshow(w_unit="nm", ax=ax[0])
    SimulatedRainbow(dw=0.2 * u.micron).imshow(ax=ax[1], w_unit="nm")
    plt.savefig(os.path.join(test_directory, "imshow-demonstration.pdf"))


def test_plot():
    SimulatedRainbow(R=10).plot()
    plt.savefig(os.path.join(test_directory, "plot-demonstration.pdf"))


def test_animate():
    # test a transit, since along both dimensions
    d = SimulatedRainbow(dw=0.1 * u.micron, dt=5 * u.minute, signal_to_noise=1000)
    theta = np.linspace(0, 2 * np.pi, d.nwave)
    planet_radius = np.sin(theta) * 0.05 + 0.15
    e = d.inject_transit(planet_radius=planet_radius)
    plotkw = dict(color="black", marker="o", linewidth=0)
    e.animate_lightcurves(
        filename=os.path.join(test_directory, "animate-lightcurves-demonstration.gif"),
        plotkw=plotkw,
    )
    e.animate_spectra(
        filename=os.path.join(test_directory, "animate-spectra-demonstration.gif"),
        plotkw=plotkw,
    )
