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
