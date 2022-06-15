from ..rainbows import *
from .setup_tests import *
from ..resampling import *


def test_bin_in_time():

    s = SimulatedRainbow(dt=5 * u.minute).inject_noise(signal_to_noise=100)
    b = s.bin_in_time(dt=1 * u.hour)
    assert b.ntime < s.ntime

    assert s.bin_in_time(dt=None, time=None, time_edges=None, ntimes=None) == s

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict(vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    b.imshow(ax=ax[1], **imshowkw)
    plt.savefig(os.path.join(test_directory, "imshow-bin-time-demonstration.pdf"))


def test_bin_in_wavelength():

    s = SimulatedRainbow(dw=50 * u.nm).inject_noise()
    b = s.bin_in_wavelength(dw=500 * u.nm)
    assert b.nwave < s.nwave

    assert (
        s.bin_in_wavelength(
            R=None, dw=None, wavelength=None, wavelength_edges=None, nwavelengths=None
        )
        == s
    )

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict(vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    b.imshow(ax=ax[1], **imshowkw)
    plt.savefig(os.path.join(test_directory, "imshow-bin-wavelength-demonstration.pdf"))
    plt.close("all")


def test_bin():

    s = SimulatedRainbow(dt=10 * u.minute, R=25).inject_noise()
    b = s.bin(dt=0.5 * u.hour, dw=0.5 * u.micron)
    assert b.nwave < s.nwave
    assert b.ntime < s.ntime

    assert (
        s.bin(
            dw=None,
            dt=None,
            R=None,
            time=None,
            time_edges=None,
            ntimes=None,
            wavelength=None,
            wavelength_edges=None,
            nwavelengths=None,
        )
        == s
    )

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict(vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    plt.title("Unbinned")
    b.imshow(ax=ax[1], **imshowkw)
    plt.title("Binned")
    plt.tight_layout()
    plt.savefig(os.path.join(test_directory, "imshow-bin-demonstration.pdf"))
    plt.close("all")


def test_bin_input_options():
    r = SimulatedRainbow(dw=0.2 * u.micron).inject_noise()

    a = r.bin(nwavelengths=3)
    centers = a.wavelength
    b = r.bin(wavelength=centers)
    edges = leftright_to_edges(a.wavelength_lower, a.wavelength_upper)
    c = r.bin(wavelength_edges=edges)
    dw = a.wavelength_upper - a.wavelength_lower
    for x in [a, b, c]:
        print(x.wavelength_lower)
        for k in r.wavelike:
            assert np.all(np.isclose(x.wavelike[k], a.wavelike[k]))

    a = r.bin(dw=0.47 * u.micron)
    centers = a.wavelength
    b = r.bin(wavelength=centers)
    edges = leftright_to_edges(a.wavelength_lower, a.wavelength_upper)
    c = r.bin(wavelength_edges=edges)
    dw = a.wavelength_upper - a.wavelength_lower
    for x in [a, b, c]:
        print(x.wavelength_lower)
        for k in r.wavelike:
            assert np.all(np.isclose(x.wavelike[k], a.wavelike[k]))

    a = r.bin(ntimes=5)
    centers = a.time
    b = r.bin(time=centers)
    edges = leftright_to_edges(a.time_lower, a.time_upper)
    c = r.bin(time_edges=edges)
    dw = a.time_upper - a.time_lower
    for x in [a, b, c]:
        print(x.time_lower)
        for k in r.timelike:
            assert np.all(np.isclose(x.timelike[k], a.timelike[k]))

    a = r.bin(dt=0.37 * u.hour)
    centers = a.time
    b = r.bin(time=centers)
    edges = leftright_to_edges(a.time_lower, a.time_upper)
    c = r.bin(time_edges=edges)
    dw = a.time_upper - a.time_lower
    for x in [a, b, c]:
        print(x.time_lower)
        for k in r.timelike:
            assert np.all(np.isclose(x.timelike[k], a.timelike[k]))


def test_bin_uncertainty_basic(original_resolution=100, binned_resolution=42):

    a = SimulatedRainbow(R=original_resolution).inject_noise()
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
    s = SimulatedRainbow(R=20, dt=5 * u.minute).inject_noise(signal_to_noise=100)
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
    plt.close("all")
    return s


def test_bin_both():
    bintime = 5
    binwave = 0.1
    w = np.logspace(0, 1) * u.micron
    r = SimulatedRainbow(dt=bintime * u.minute).inject_noise(signal_to_noise=1000)
    b_withouttransit = r.bin(dw=binwave * u.micron, dt=bintime * u.minute)


def test_binning_to_one():
    """
    Are things OK even if we bin down to one wavelength or time?
    """
    for N in np.random.randint(2, 100, 3):
        for M in np.random.randint(2, 100, 3):
            w = np.linspace(1, 2, N) * u.micron
            t = np.linspace(-1, 1, M) * u.hour
            s = SimulatedRainbow(wavelength=w, time=t).inject_noise()

            bw = s.bin(nwavelengths=s.nwave)
            print(bw.uncertainty[0, 0], s.uncertainty[0, 0] / np.sqrt(s.nwave))
            predicted = 1 / np.sqrt(np.sum(1 / s.uncertainty**2, axis=s.waveaxis))
            assert np.all(
                np.isclose(
                    predicted,
                    np.mean(bw.uncertainty, axis=s.waveaxis),
                    rtol=1e-10,
                    atol=1e-13,
                )
            )
            assert bw.nwave == 1

            bt = s.bin(ntimes=s.ntime)
            print(bt.uncertainty[0, 0], s.uncertainty[0, 0] / np.sqrt(s.ntime))
            predicted = 1 / np.sqrt(np.sum(1 / s.uncertainty**2, axis=s.timeaxis))
            assert np.all(
                np.isclose(
                    predicted,
                    np.mean(bt.uncertainty, axis=s.timeaxis),
                    rtol=1e-10,
                    atol=1e-13,
                )
            )
            assert bt.ntime == 1


def test_integrated_wrappers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        s = SimulatedRainbow().inject_transit().inject_noise()
        s.fluxlike["flux"] += 0.003 * s.wavelength.value[:, np.newaxis]

        fi, ax = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
        kw = dict(vmin=0.98, vmax=1.02)
        s.imshow(ax=ax[0], **kw)
        s.get_spectrum_as_rainbow().imshow(ax=ax[1], **kw)
        s.get_lightcurve_as_rainbow().imshow(ax=ax[2], **kw)
