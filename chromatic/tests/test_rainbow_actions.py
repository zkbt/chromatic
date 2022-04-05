from ..rainbows import *
from .setup_tests import *


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


def test_bin_uncertainty_basic(original_resolution=100, binned_resolution=42):

    a = SimulatedRainbow(R=original_resolution, signal_to_noise=100)
    b = a.bin(R=binned_resolution)
    original_uncertainty = np.median(a.uncertainty)
    predicted_uncertainty = np.median(a.uncertainty) / np.sqrt(
        original_resolution / binned_resolution
    )
    actual_uncertainty = np.median(b.uncertainty)

    did_it_work = np.isclose(predicted_uncertainty, actual_uncertainty, rtol=0.001)

    print(
        f"""
    At the original resolution of R={original_resolution},
    the median uncertainty was {original_uncertainty}.
    When binned to R={binned_resolution}, the predicted
    resolution should be {predicted_uncertainty}.
    It is actually {actual_uncertainty},
    so it is {did_it_work} that the binning worked!
    """
    )

    assert did_it_work


def test_bin_bad_data(visualize=False):
    s = SimulatedRainbow(signal_to_noise=100, R=20, dt=5 * u.minute)
    N = 100
    i, j = np.random.randint(0, s.nwave, N), np.random.randint(0, s.ntime, N)
    s.fluxlike["flux"][i, j] = np.nan

    i, j = np.random.randint(0, s.nwave, N), np.random.randint(0, s.ntime, N)
    s.fluxlike["uncertainty"][i, j] = 0

    b = s.bin(R=7, dt=15 * u.minute)

    if visualize:
        fi, ax = plt.subplots(2, 2, figsize=(8, 6), dpi=300, sharex=True, sharey=True)
        s.imshow(ax=ax[0, 0], quantity="uncertainty", vmax=0.01)
        s.imshow(ax=ax[0, 1], vmin=0.97, vmax=1.03)

        b.imshow(ax=ax[1, 0], quantity="uncertainty", vmax=0.01)
        b.imshow(ax=ax[1, 1], vmin=0.97, vmax=1.03)

    assert np.any(np.isfinite(b.flux))
    return s


def test_normalize(plot=False):
    N = 37
    w = np.logspace(0, 1, N) * u.micron
    f = np.cos(2 * np.pi * w.value / 3) + 1

    snr = 100
    a = SimulatedRainbow(signal_to_noise=snr, wavelength=w)
    b = SimulatedRainbow(signal_to_noise=snr, wavelength=w, star_flux=f)
    c = SimulatedRainbow(
        signal_to_noise=snr, wavelength=w, star_flux=f
    ).inject_transit()

    for x in [a, b, c]:
        nw = x.normalize(axis="w")
        nt = x.normalize(axis="t")
        nwt = x.normalize(axis="w").normalize(axis="t")
        ntw = x.normalize(axis="t").normalize(axis="w")

        for r in [nw, nt, nwt, ntw]:
            r.fluxlike["relative-uncertainty"] = r.uncertainty / r.flux
            assert np.all(np.isclose(r.uncertainty / r.flux, 1 / snr, rtol=0.1))
            if plot:
                r.imshow_fluxlike_quantities(maxcol=4)
