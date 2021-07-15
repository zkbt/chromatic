from ..rainbows import *
from .setup_tests import *


def test_simulated_basics():

    a = SimulatedRainbow()
    b = SimulatedRainbow(tlim=[-1, 1] * u.hour, dt=1 * u.minute)
    c = SimulatedRainbow(time=np.linspace(-1, 1) * u.hour)
    d = SimulatedRainbow(wlim=[0.1, 1] * u.micron, R=100)
    assert d.wscale == "log"
    e = SimulatedRainbow(wlim=[0.1, 1] * u.micron, dw=1 * u.nm)
    assert e.wscale == "linear"
    f = SimulatedRainbow(wavelength=np.logspace(0, 1) * u.micron)


def test_star_flux():
    g = SimulatedRainbow(
        wavelength=np.logspace(0, 1) * u.micron, star_flux=np.logspace(0, 1)
    )


def test_inject_transit():
    h = SimulatedRainbow(wavelength=np.logspace(0, 1) * u.micron)
    h.inject_transit(planet_params={"per": 0.1})
    h.inject_transit(planet_radius=np.zeros(50) + 0.1)

    i = SimulatedRainbow(signal_to_noise=1000, dt=2 * u.minute)
    fi, ax = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    i.imshow(ax=ax[0], vmin=0.975, vmax=1.005)
    plt.xlabel("")
    i.inject_transit(
        planet_params=dict(per=3), planet_radius=np.random.normal(0.1, 0.01, i.nwave)
    ).imshow(ax=ax[1], vmin=0.975, vmax=1.005)
    # test the limb darkening
    i.inject_transit(
        planet_params={
            "per": 3,
            "limb_dark": "quadratic",
            "u": np.transpose(
                [np.linspace(1.0, 0.0, i.nwave), np.linspace(0.5, 0.0, i.nwave)]
            ),
        },
    ).imshow(ax=ax[2], vmin=0.975, vmax=1.005)
    plt.savefig(os.path.join(test_directory, "transit-injection-demonstration.pdf"))
