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
    scatterkw = dict()
    e.animate_lightcurves(
        filename=os.path.join(test_directory, "animate-lightcurves-demonstration.gif"),
        scatterkw=scatterkw,
    )
    e.animate_spectra(
        filename=os.path.join(test_directory, "animate-spectra-demonstration.gif"),
        scatterkw=scatterkw,
    )


def test_wavelength_cmap():

    r = SimulatedRainbow(R=10)

    # can we set up the wavelength-based color map
    r._setup_wavelength_colors(cmap=one2another("black", "red"))

    # test a few examples
    assert r.get_wavelength_color(r.wavelength[0]) == (0.0, 0.0, 0.0, 1.0)
    assert r.get_wavelength_color(r.wavelength[-1]) == (1.0, 0.0, 0.0, 1.0)
