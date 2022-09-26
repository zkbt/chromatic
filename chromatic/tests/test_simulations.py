from ..rainbows import *
from .setup_tests import *


def test_simulated_basics():

    a = SimulatedRainbow().inject_noise()
    b = SimulatedRainbow(tlim=[-1, 1] * u.hour, dt=1 * u.minute).inject_noise()
    c = SimulatedRainbow(time=np.linspace(-1, 1) * u.hour).inject_noise()
    d = SimulatedRainbow(wlim=[0.1, 1] * u.micron, R=100).inject_noise()
    assert d.wscale == "log"
    e = SimulatedRainbow(wlim=[0.1, 1] * u.micron, dw=1 * u.nm).inject_noise()
    assert e.wscale == "linear"
    f = SimulatedRainbow(wavelength=np.logspace(0, 1) * u.micron).inject_noise()


def test_photon_noise(N=10000):

    # inject noise in two different ways
    noiseless = SimulatedRainbow().inject_transit()
    gaussian = noiseless.inject_noise(signal_to_noise=np.sqrt(N)).normalize()
    poisson = noiseless.inject_noise(number_of_photons=N).normalize()

    # imshow the two ways
    fi, ax = plt.subplots(1, 2, figsize=(8, 3), dpi=300, constrained_layout=True)
    gaussian.imshow(ax=ax[0])
    ax[0].set_title(r"Gaussian ($\sigma=$" + f"{1/np.sqrt(N):.3f}=" + r"$1/\sqrt{N}$)")
    poisson.imshow(ax=ax[1])
    ax[1].set_title(f"Poisson (N={N})")

    # make sure the standard deviation is about right
    # sigma = np.std(poisson.residuals / poisson.model)
    # assert np.isclose(sigma, 1 / np.sqrt(N), rtol=0.1)
    plt.savefig(
        os.path.join(
            test_directory,
            "demonstration-of-injecting-noise-as-poisson-or-gaussian.pdf",
        )
    )


def test_star_flux():
    g = SimulatedRainbow(
        wavelength=np.logspace(0, 1) * u.micron, star_flux=np.logspace(0, 1)
    ).inject_noise()


def test_inject_transit():
    h = SimulatedRainbow(wavelength=np.logspace(0, 1) * u.micron).inject_noise()
    h.inject_transit(per=0.1)
    h.inject_transit(planet_radius=np.zeros(50) + 0.1)

    i = SimulatedRainbow(dt=2 * u.minute).inject_noise(signal_to_noise=1000)
    fi, ax = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    i.imshow(ax=ax[0], vmin=0.975, vmax=1.005)
    plt.xlabel("")
    i.inject_transit(per=3, planet_radius=np.random.normal(0.1, 0.01, i.nwave)).imshow(
        ax=ax[1], vmin=0.975, vmax=1.005
    )
    # test the limb darkening
    i.inject_transit(
        **{
            "per": 3,
            "limb_dark": "quadratic",
            "u": np.transpose(
                [np.linspace(1.0, 0.0, i.nwave), np.linspace(0.5, 0.0, i.nwave)]
            ),
        },
    ).imshow(ax=ax[2], vmin=0.975, vmax=1.005)
    plt.savefig(os.path.join(test_directory, "demonstration-of-injecting-transit.pdf"))
    plt.close("all")


def test_inject_systematics():
    SimulatedRainbow().inject_transit().inject_noise().inject_systematics().bin(
        R=10, dt=10 * u.minute
    ).imshow_with_models()
    plt.savefig(
        os.path.join(test_directory, "demonstration-of-injecting-fake-systematics.pdf")
    )
