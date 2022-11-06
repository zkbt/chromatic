from ..rainbows import *
from .setup_tests import *


def test_attach_model():
    plt.close("all")

    simulated = SimulatedRainbow().inject_noise().inject_transit().inject_systematics()
    original = simulated._create_copy()

    # pull out the model
    inputs = {}
    for k in ["model", "planet_model", "systematics_model"]:
        inputs[k] = simulated.fluxlike[k]
        del simulated.fluxlike[k]

    # reattach the model
    withmodel = simulated.attach_model(**inputs)
    assert simulated != withmodel
    assert withmodel == original


def test_imshow_with_models():
    plt.close("all")

    s = (
        SimulatedRainbow()
        .inject_transit()
        .inject_systematics()
        .inject_noise()
        .bin(R=50, dt=5 * u.minute)
    )
    s.imshow_with_models(cmap="gray")
    plt.savefig(
        os.path.join(test_directory, "demonstration-of-imshow-data-with-model.pdf")
    )
    s.imshow_with_models(models=["systematics_model", "planet_model"], cmap="gray")
    plt.savefig(
        os.path.join(
            test_directory, "demonstration-of-imshow-data-with-model-components.pdf"
        )
    )

    s.imshow_with_models(models=["systematics_model", "planet_model"], cmap="gray")
    s.imshow_with_models(
        models=["systematics_model", "planet_model"], cmap="gray", label=False
    )
    s.imshow_with_models(
        models=["systematics_model", "planet_model"],
        cmap="gray",
        label_textkw=dict(color="red"),
    )
    s.imshow_with_models(
        models=["systematics_model", "planet_model"], cmap="gray", label="outside"
    )


def test_plot_with_model_and_residuals():
    plt.close("all")

    s = SimulatedRainbow(R=3, dt=3 * u.minute)
    r = s.inject_transit(
        limb_dark="quadratic",
        u=np.transpose(
            [np.linspace(1.0, 0.0, s.nwave), np.linspace(0.5, 0.0, s.nwave)]
        ),
        method="batman",
    ).inject_noise(signal_to_noise=1000)
    for i, options in enumerate(
        [
            dict(),
            dict(errorbar=True),
            dict(cmap=one2another("skyblue", "sienna")),
            dict(data_plotkw=dict(alpha=0.5), cmap="magma_r"),
            dict(figsize=(8, 4), cmap=one2another("orchid", "indigo")),
        ]
    ):
        r.plot_with_model_and_residuals(**options)
        plt.savefig(
            os.path.join(
                test_directory,
                f"demonstration-of-rainbow-of-lightcurves-and-residuals-example{i}.png",
            )
        )


def test_plot_and_animate_with_models(output="gif"):
    plt.close("all")
    s = (
        SimulatedRainbow(R=5)
        .inject_transit()
        .inject_systematics()
        .inject_noise(signal_to_noise=500)
    )
    s.setup_wavelength_colors(
        cmap="magma", vmin=0.1 * u.micron, vmax=10 * u.micron, log=True
    )
    for o in ["vertical", "horizontal"]:
        for e in [True, False]:
            s.plot_one_wavelength_with_models(0, errorbar=e, orientation=o)
            error_string = {True: "with", False: "without"}[e] + "-errorbars"
            filename = f"data-with-models-{o}-{error_string}"
            plt.savefig(
                os.path.join(test_directory, f"demonstration-of-plot-{filename}.png")
            )
            s.animate_with_models(
                os.path.join(
                    test_directory, f"demonstration-of-animate-{filename}.{output}"
                ),
                errorbar=e,
                orientation=o,
            )
