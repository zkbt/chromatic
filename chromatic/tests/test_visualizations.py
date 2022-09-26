from ..rainbows import *
from .setup_tests import *
from ..rainbows.visualizations import _add_panel_labels


def test_imshow():
    plt.close("all")

    plt.figure()
    SimulatedRainbow(R=20).imshow()
    plt.savefig(os.path.join(test_directory, "demonstration-of-imshow.pdf"))

    plt.figure()
    SimulatedRainbow(dw=0.2 * u.micron).imshow()

    plt.figure()
    Rainbow(
        wavelength=np.arange(1, 5) * u.micron,
        time=np.arange(1, 6) * u.hour,
        flux=np.ones((4, 5)),
    ).imshow()

    fi, ax = plt.subplots(2, 1, sharex=True)
    SimulatedRainbow(R=10).inject_noise().imshow(w_unit="nm", ax=ax[0])
    SimulatedRainbow(dw=0.2 * u.micron).inject_noise().imshow(ax=ax[1], w_unit="nm")
    plt.savefig(os.path.join(test_directory, "demonstration-of-imshow-units.pdf"))


def test_imshow_quantities():
    plt.close("all")

    s = SimulatedRainbow().inject_transit().inject_noise(signal_to_noise=500)
    for k in "abcde":
        s.fluxlike[k] = np.random.uniform(4, 5, s.shape)
    s.imshow_quantities(maxcol=1, panel_size=(8, 2))
    plt.savefig(os.path.join(test_directory, "demonstration-of-imshow-quantities.pdf"))


def test_plot():
    plt.close("all")

    SimulatedRainbow(R=10).inject_noise().plot()
    plt.savefig(os.path.join(test_directory, "demonstration-of-plot.pdf"))


def test_plot_unnormalized():
    plt.close("all")

    w = np.logspace(0, 1, 5) * u.micron
    plt.figure()
    s = SimulatedRainbow(wavelength=w, star_flux=w.value**2).inject_noise(
        signal_to_noise=5
    )
    s.plot(spacing=0)
    plt.savefig(
        os.path.join(
            test_directory,
            "demonstration-of-plot-without-normalization.pdf",
        )
    )


def test_plot_quantities():
    plt.close("all")

    r = SimulatedRainbow(R=10).inject_noise()
    for k in "abcdefg":
        r.timelike[f'timelike quantity "{k}"'] = np.random.normal(0, 1, r.ntime) * u.m
        r.wavelike[f'wavelike quantity "{k}"'] = np.random.normal(0, 1, r.nwave) * u.s

    for k in ["time", "wavelength"]:
        for x in [k, "index"]:
            r.plot_quantities(xaxis=k, what_is_x=x)
            plt.savefig(
                os.path.join(
                    test_directory,
                    f"demonstration-of-plot-quantities-data={k}-xaxis={x}.pdf",
                )
            )


def test_animate():
    plt.close("all")

    # test a transit, since along both dimensions
    d = SimulatedRainbow(dw=0.1 * u.micron, dt=5 * u.minute).inject_noise(
        signal_to_noise=1000
    )
    theta = np.linspace(0, 2 * np.pi, d.nwave)
    planet_radius = np.sin(theta) * 0.05 + 0.15
    e = d.inject_transit(planet_radius=planet_radius)
    scatterkw = dict()
    e.animate_lightcurves(
        filename=os.path.join(
            test_directory, "demonstration-of-animate-lightcurves.gif"
        ),
        scatterkw=scatterkw,
    )
    e.animate_spectra(
        filename=os.path.join(
            test_directory, "demonstration-of-animate-spectra-demonstration.gif"
        ),
        scatterkw=scatterkw,
    )


def test_animate_other_quantities():
    plt.close("all")

    k = "some-imaginary-fluxlike-quantity"
    s = (
        SimulatedRainbow(R=5, dt=20 * u.minute)
        .inject_transit()
        .inject_systematics(amplitude=0.001, fluxlike=[k])
        .inject_noise(signal_to_noise=1000)
    )
    s.animate_spectra(
        os.path.join(
            test_directory,
            "demonstration-of-animate-quantities-beside-flux-spectra.gif",
        ),
        quantity=k,
    )
    s.animate_lightcurves(
        os.path.join(
            test_directory,
            "demonstration-of--animate-quantities-beside-flux-lightcurves.gif",
        ),
        quantity=k,
    )


def test_cmap():
    plt.close("all")

    r = SimulatedRainbow(R=10).inject_noise()

    # can we set up the wavelength-based color map
    r.setup_wavelength_colors(cmap=one2another("black", "red"))

    # test a few examples
    assert r.get_wavelength_color(r.wavelength[0]) == (0.0, 0.0, 0.0, 1.0)
    assert r.get_wavelength_color(r.wavelength[-1]) == (1.0, 0.0, 0.0, 1.0)


def test_imshow_interact():
    plt.close("all")

    plt.figure()
    SimulatedRainbow(R=10).inject_noise().imshow_interact()


def test_plot_one_wavelength():
    plt.close("all")

    s = SimulatedRainbow(wavelength=[1] * u.micron).inject_noise()
    s.plot()


def test_imshow_one_wavelength():
    plt.close("all")

    with pytest.warns(match="hard to imshow "):

        s = SimulatedRainbow(wavelength=[1] * u.micron).inject_noise()
        ax = s.imshow()
        assert "Wavelength Index" in ax.get_ylabel()

        s = SimulatedRainbow().inject_noise()
        b = s.bin(nwavelengths=s.nwave)
        ax = b.imshow()
        assert "Wavelength (" in ax.get_ylabel()
        ylim = ax.get_ylim()
    plt.savefig(
        os.path.join(
            test_directory,
            f"demonstration-of-imshow-with-just-one-wavelength.pdf",
        )
    )


def test_imshow_randomized_axes():
    plt.close("all")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        s = (
            SimulatedRainbow(
                time=np.random.uniform(-3, 3, 12) * u.hour,
                wavelength=np.random.uniform(0.5, 5, 8) * u.micron,
            )
            .inject_transit()
            .inject_noise(signal_to_noise=300)
        )
        s.fluxlike["flux"] += 0.003 * s.wavelength.value[:, np.newaxis]

        fi, ax = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
        kw = dict(vmin=0.98, vmax=1.02)
        s.imshow(ax=ax[0], **kw)
        s.get_average_spectrum_as_rainbow().imshow(ax=ax[1], **kw)
        s.get_average_lightcurve_as_rainbow().imshow(ax=ax[2], **kw)
        for i in [0, 2]:
            assert "Time Index" in ax[i].get_xlabel()
        for i in [0, 1]:
            assert "Wavelength Index" in ax[1].get_ylabel()
        plt.savefig(
            os.path.join(
                test_directory,
                f"demonstration-of-imshow-with-irregular-axes.pdf",
            )
        )


def test_imshow_both_orientations():
    plt.close("all")

    s = (
        SimulatedRainbow(R=5, dt=10 * u.minute)
        .inject_transit(planet_radius=0.2)
        .inject_systematics(amplitude=0.005)
        .inject_noise(signal_to_noise=1000)
    )

    fi, ax = plt.subplots(2, 2, constrained_layout=True, figsize=(6, 6))
    for i, k in enumerate(["time", "wavelength"]):
        s.plot(ax=ax[0, i], xaxis=k)
        s.imshow(ax=ax[1, i], xaxis=k, colorbar=False)
        plt.savefig(
            os.path.join(
                test_directory,
                f"demonstration-of-imshow-with-both-orientations.pdf",
            )
        )


def test_pcolormesh_both_orientations():
    plt.close("all")

    s = (
        SimulatedRainbow(R=5, dt=10 * u.minute)
        .inject_transit(planet_radius=0.2)
        .inject_systematics(amplitude=0.005)
        .inject_noise(signal_to_noise=1000)
    )

    fi, ax = plt.subplots(2, 2, constrained_layout=True, figsize=(6, 6))
    for i, k in enumerate(["time", "wavelength"]):
        s.plot(ax=ax[0, i], xaxis=k)
        s.pcolormesh(ax=ax[1, i], xaxis=k, colorbar=False)
    plt.savefig(
        os.path.join(
            test_directory, "demonstration-of-pcolormesh-both-orientations.pdf"
        )
    )


def test_both_types_of_plot():
    plt.close("all")

    N, M = 10, 20
    r = (
        SimulatedRainbow(
            wavelength=np.linspace(1, 2, N) * u.micron,
            time=np.linspace(-3, 3, M) * u.hour,
            star_flux=np.linspace(1, 1.1, N),
        )
        .inject_transit(planet_radius=0.2)
        .inject_noise(signal_to_noise=500)
    )
    fi, ax = plt.subplots(4, 2, figsize=(8, 15), constrained_layout=True)
    for i, f in enumerate(["plot_lightcurves", "plot_spectra"]):
        getattr(r, f)(ax=ax[0, i])
        getattr(r.normalize(), f)(ax=ax[1, i], errorbar=True, spacing=0)
        getattr(r.normalize(), f)(ax=ax[2, i], spacing=0.01)
        getattr(r.normalize(), f)(
            ax=ax[3, i],
            spacing=0,
            errorbar=True,
            plotkw=dict(color="orchid"),
            textkw=dict(color="orchid"),
            errorbarkw=dict(color="orchid"),
        )
    plt.savefig(
        os.path.join(
            test_directory, "demonstration-of-plot-lightcurves-and-plot-spectra.pdf"
        )
    )


def test_add_labels_to_panels():
    plt.close("all")

    fi, ax = plt.subplots(3, 3)
    _add_panel_labels(ax, preset="inside", color="blue")
    plt.savefig(
        os.path.join(
            test_directory, "demonstration-of-adding-abcde-labels-to-panels.pdf"
        )
    )


def test_pcolormesh():
    plt.close("all")

    s = (
        SimulatedRainbow(R=10, dt=10 * u.minute)
        .inject_transit()
        .inject_noise(signal_to_noise=1000)
    )
    fi, ax = plt.subplots(2, 2, constrained_layout=True, figsize=(8, 5))
    s.imshow(ax=ax[0, 0])
    plt.title("imshow")
    s.pcolormesh(ax=ax[0, 1])
    plt.title("pcolormesh")
    plt.yscale("log")
    b = s.bin(dw=0.5 * u.micron)
    b.imshow(ax=ax[1, 0])
    b.pcolormesh(ax=ax[1, 1])
    plt.savefig(
        os.path.join(test_directory, "demonstration-of-pcolormesh-vs-imshow.pdf")
    )


def test_plot_noise_comparison():
    fi, ax = plt.subplots(
        2,
        2,
        sharey="row",
        figsize=(8, 6),
        sharex=True,
        dpi=300,
        constrained_layout=True,
    )
    for i, a in enumerate([0, 0.01]):
        s = (
            SimulatedRainbow(dw=0.1 * u.micron)
            .inject_systematics(amplitude=a)
            .inject_noise()
        )
        s.imshow(ax=ax[0, i], xaxis="wavelength")
        s.plot_noise_comparison(ax=ax[1, i])
    ax[0, 0].set_title("No Systematics")
    ax[0, 1].set_title("With Systematics")
    plt.savefig(
        os.path.join(test_directory, "demonstration-of-plot-noise-comparison.pdf")
    )


def test_plot_noise_comparison_in_bins():
    plt.close("all")

    s = SimulatedRainbow(dw=0.1 * u.micron).inject_systematics().inject_noise()
    s.plot_noise_comparison_in_bins()
    plt.savefig(
        os.path.join(
            test_directory, "demonstration-of-plot-noise-comparison-in-bins.pdf"
        )
    )


def test_plot_histogram():
    plt.close("all")

    s = SimulatedRainbow(R=5).inject_noise()
    fi, ax = plt.subplots(
        s.nwave, 1, figsize=(4, 12), sharex=True, sharey=True, constrained_layout=True
    )
    for i in range(s.nwave):
        s.plot_histogram(i, ax=ax[i], expected=True)
        if i < (s.nwave - 1):
            plt.xlabel("")
    plt.savefig(os.path.join(test_directory, "demonstration-of-histogram.pdf"))
