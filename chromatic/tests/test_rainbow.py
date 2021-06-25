from ..rainbows import *
from .setup_tests import *
import pytest


def test_basic_rainbow():

    nw, nt = 20, 40
    r = Rainbow(
        wavelength=np.linspace(0.5, 5, nw) * u.micron,
        time=np.linspace(-1, 1, nt) * u.hour,
        flux=np.ones((nw, nt)),
    )

    assert r.shape == (nw, nt)
    assert r.nwave == nw
    assert r.ntime == nt
    assert r.nflux == nw * nt


def test_bin_in_time():

    s = SimulatedRainbow(dt=5 * u.minute, signal_to_noise=100)
    b = s.bin_in_time(dt=1 * u.hour)
    assert b.ntime < s.ntime

    assert s.bin_in_time(dt=None, time=None) == s

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict(vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    b.imshow(ax=ax[1], **imshowkw)
    plt.savefig(os.path.join(test_directory, "imshow-bin-time-demonstration.pdf"))


def test_bin_in_wavelength():

    s = SimulatedRainbow(dw=50 * u.nm)
    b = s.bin_in_wavelength(dw=500 * u.nm)
    assert b.nwave < s.nwave

    assert s.bin_in_wavelength(R=None, dw=None, wavelength=None) == s

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict(vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    b.imshow(ax=ax[1], **imshowkw)
    plt.savefig(os.path.join(test_directory, "imshow-bin-wavelength-demonstration.pdf"))


def test_bin():

    s = SimulatedRainbow(dt=10 * u.minute, R=25)
    b = s.bin(dt=0.5 * u.hour, dw=0.5 * u.micron)
    assert b.nwave < s.nwave
    assert b.ntime < s.ntime

    assert s.bin(dw=None, dt=None, R=None, time=None, wavelength=None) == s

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict(vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    plt.title("Unbinned")
    b.imshow(ax=ax[1], **imshowkw)
    plt.title("Binned")
    plt.tight_layout()
    plt.savefig(os.path.join(test_directory, "imshow-bin-demonstration.pdf"))


def test_rainbow_operations():
    nw, nt = 20, 40
    a = Rainbow(
        wavelength=np.linspace(0.5, 5, nw) * u.micron,
        time=np.linspace(-1, 1, nt) * u.hour,
        flux=np.ones((nw, nt)),
    )

    b = Rainbow(
        wavelength=np.linspace(0.5, 5, nw) * u.micron,
        time=np.linspace(-1, 1, nt) * u.hour,
        flux=np.zeros((nw, nt)) + 0.1,
    )

    wl_like = np.linspace(-0.01, 0.01, nw)
    t_like = np.linspace(-0.01, 0.01, nt)
    fx_like = np.ones(a.shape)

    assert (a + b).fluxlike["flux"][0][0] == 1.1
    assert (a + wl_like).fluxlike["flux"][0][0] == 0.99
    assert (a + t_like).fluxlike["flux"][0][0] == 0.99
    assert (a + fx_like).fluxlike["flux"][0][0] == 2
    assert (a + 1).fluxlike["flux"][0][0] == 2

    assert (a - b).fluxlike["flux"][0][0] == 0.9
    assert (a - wl_like).fluxlike["flux"][0][0] == 1.01
    assert (a - t_like).fluxlike["flux"][0][0] == 1.01
    assert (a - fx_like).fluxlike["flux"][0][0] == 0
    assert (a - 1).fluxlike["flux"][0][0] == 0

    assert (a * b).fluxlike["flux"][0][0] == 0.1
    assert (a * wl_like).fluxlike["flux"][0][0] == -0.01
    assert (a * t_like).fluxlike["flux"][0][0] == -0.01
    assert (a * fx_like).fluxlike["flux"][0][0] == 1
    assert (a * 1).fluxlike["flux"][0][0] == 1

    assert (a / b).fluxlike["flux"][0][0] == 1 / 0.1
    assert (a / wl_like).fluxlike["flux"][0][0] == 1 / -0.01
    assert (a / t_like).fluxlike["flux"][0][0] == 1 / -0.01
    assert (a / fx_like).fluxlike["flux"][0][0] == 1
    assert (a / 1).fluxlike["flux"][0][0] == 1

    # make sure we raise an error if it's not obvious whether we're doing wavelength or time
    c = Rainbow(
        wavelength=np.linspace(0.5, 5, nw) * u.micron,
        time=np.linspace(-1, 1, nw) * u.hour,
        flux=np.ones((nw, nw)),
    )
    with pytest.raises(Exception):
        c * wl_like
    with pytest.raises(Exception):
        c / wl_like
    with pytest.raises(Exception):
        c + wl_like
    with pytest.raises(Exception):
        c - wl_like


def test_animate():
    star_wl = np.arange(0.5, 5, 0.2)
    star_fx = np.ones(len(star_wl)) + np.linspace(-0.01, 0.01, len(star_wl))
    planet_radius = np.linspace(0.1, 0.2, len(star_wl))
    d = SimulatedRainbow(
        dw=0.2 * u.micron, star_flux=star_fx, dt=10 * u.minute, signal_to_noise=1000
    )
    e = d.inject_transit(planet_radius=planet_radius)
    pltkw = dict(color="black", marker="o", linewidth=0)
    e.animate_lightcurve(
        os.path.join(test_directory, "animate-lightcurve-demonstration.gif"), **pltkw
    )
    e.animate_spectra(
        os.path.join(test_directory, "animate-spectra-demonstration.gif"), **pltkw
    )
