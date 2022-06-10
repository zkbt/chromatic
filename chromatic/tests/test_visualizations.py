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


def test_imshow_quantities():
    s = SimulatedRainbow(signal_to_noise=500).inject_transit()
    for k in "abcde":
        s.fluxlike[k] = np.random.uniform(4, 5, s.shape)
    s.imshow_quantities(maxcol=1, panel_size=(8, 2))
    plt.savefig(os.path.join(test_directory, "imshow-multiples-demonstration.pdf"))


def test_plot():
    SimulatedRainbow(R=10).plot()
    plt.savefig(os.path.join(test_directory, "plot-demonstration.pdf"))


def test_plot_unnormalized():
    w = np.logspace(0, 1, 5) * u.micron
    plt.figure()
    s = SimulatedRainbow(wavelength=w, star_flux=w.value**2, signal_to_noise=5)
    s.plot(spacing=0)
    plt.savefig(os.path.join(test_directory, "plot-demonstration-unnormalized.pdf"))


def test_plot_quantities():
    r = SimulatedRainbow(R=10)
    for k in "abcdefg":
        r.timelike[f'timelike quantity "{k}"'] = np.random.normal(0, 1, r.ntime) * u.m
        r.wavelike[f'wavelike quantity "{k}"'] = np.random.normal(0, 1, r.nwave) * u.s

    for k in ["time", "wavelength"]:
        for x in [k, "index"]:
            r.plot_quantities(data_like=k, x_axis=x)
            plt.savefig(
                os.path.join(
                    test_directory,
                    f"plot_quantities-demonstration-data={k}-xaxis={x}.pdf",
                )
            )


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


def test_imshow_interact():
    plt.figure()
    SimulatedRainbow(R=10).imshow_interact()


def test_plot_one_wavelength():
    s = SimulatedRainbow(wavelength=[1] * u.micron)
    s.plot()


def test_imshow_one_wavelength():
    s = SimulatedRainbow(wavelength=[1] * u.micron)
    ax = s.imshow()
    assert "Wavelength Index" in ax.get_ylabel()

    s = SimulatedRainbow()
    b = s.bin(nwavelengths=s.nwave)
    ax = b.imshow()
    assert "Wavelength (" in ax.get_ylabel()
    ylim = ax.get_ylim()
    assert ylim[0] > max(s.wavelength.value)
    assert ylim[1] < min(s.wavelength.value)
    plt.close("all")


def test_imshow_randomized_axes():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        s = SimulatedRainbow(
            time=np.random.uniform(-3, 3, 12) * u.hour,
            wavelength=np.random.uniform(0.5, 5, 8) * u.micron,
            signal_to_noise=300,
        ).inject_transit()
        s.fluxlike["flux"] += 0.003 * s.wavelength.value[:, np.newaxis]

        fi, ax = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
        kw = dict(vmin=0.98, vmax=1.02)
        s.imshow(ax=ax[0], **kw)
        s.get_spectrum_as_rainbow().imshow(ax=ax[1], **kw)
        s.get_lightcurve_as_rainbow().imshow(ax=ax[2], **kw)
        for i in [0, 2]:
            assert "Time Index" in ax[i].get_xlabel()
        for i in [0, 1]:
            assert "Wavelength Index" in ax[1].get_ylabel()
