from ..rainbows import *
from .setup_tests import *


def test_attach_model():
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


def test_imshow_data_with_models():
    s = (
        SimulatedRainbow()
        .inject_transit()
        .inject_systematics()
        .inject_noise()
        .bin(R=50, dt=5 * u.minute)
    )
    s.imshow_data_with_models(cmap="gray")
    plt.savefig(os.path.join(test_directory, "imshow-data-with-model.pdf"))
    s.imshow_data_with_models(models=["systematics_model", "planet_model"], cmap="gray")
    plt.savefig(os.path.join(test_directory, "imshow-data-with-model-components.pdf"))

    s.imshow_data_with_models(models=["systematics_model", "planet_model"], cmap="gray")
    s.imshow_data_with_models(
        models=["systematics_model", "planet_model"], cmap="gray", label=False
    )
    s.imshow_data_with_models(
        models=["systematics_model", "planet_model"],
        cmap="gray",
        labelkw=dict(color="red"),
    )
    s.imshow_data_with_models(
        models=["systematics_model", "planet_model"], cmap="gray", label="outside"
    )


def test_plot_lightcurves_and_residuals_with_models():
    s = SimulatedRainbow(R=3, dt=3 * u.minute)
    r = s.inject_transit(
        limb_dark="quadratic",
        u=np.transpose(
            [np.linspace(1.0, 0.0, s.nwave), np.linspace(0.5, 0.0, s.nwave)]
        ),
    ).inject_noise(signal_to_noise=1000)
    for i, options in enumerate(
        [
            dict(),
            dict(cmap=one2another("skyblue", "sienna")),
            dict(data_plotkw=dict(alpha=0.5), cmap="magma_r"),
            dict(figsize=(8, 4), cmap=one2another("orchid", "indigo")),
        ]
    ):
        r.plot_lightcurves_and_residuals_with_models(**options)
        plt.savefig(
            os.path.join(
                test_directory, f"rainbow-of-lightcurves-and-residuals-example{i}.png"
            )
        )
