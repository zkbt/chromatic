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
    s = SimulatedRainbow()
    x = (s.wavelength / max(s.wavelength)).value

    # test that injected parameters are stored appropriately
    delta = 0.5 * np.sin(x * 10)
    t0 = 0.01
    t = s.inject_transit(delta=delta, t0=t0, method="trapezoid")
    assert t.metadata["injected_transit_parameters"]["t0"] == t0
    assert np.all(t.wavelike["injected_transit_delta"] == delta)

    # test that all methods work
    fi, ax = plt.subplots(2, 2, figsize=(8, 6), constrained_layout=True, dpi=300)

    method = "trapezoid"
    s.inject_transit(method=method).imshow(ax=ax[0, 0])
    plt.title(f"{method} | default")
    s.inject_transit(
        planet_radius=0.5 * np.sin(x * 10),
        tau=0.02,
        t0=x * 0.03,
        T=x * 0.05,
        method=method,
    ).imshow(ax=ax[0, 1])
    plt.title(f"{method} | wacky")

    method = "batman"
    plt.title(f"{method} | default")
    s.inject_transit(method=method).imshow(ax=ax[1, 0])
    plt.title(f"{method} | default")
    s.inject_transit(
        planet_radius=0.5 * np.sin(x * 10),
        t0=x * 0.03,
        inc=89,
        a=1 / x * 10,
        limb_dark="nonlinear",
        u=np.random.uniform(0, 0.5, [s.nwave, 4]),
        method=method,
    ).imshow(ax=ax[1, 1])
    plt.title(f"{method} | wacky")

    plt.savefig(os.path.join(test_directory, "demonstration-of-injecting-transit.pdf"))
    plt.close("all")


def test_inject_systematics():
    SimulatedRainbow().inject_transit().inject_noise().inject_systematics().bin(
        R=10, dt=10 * u.minute
    ).imshow_with_models()
    plt.savefig(
        os.path.join(test_directory, "demonstration-of-injecting-fake-systematics.pdf")
    )
