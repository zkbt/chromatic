from ..modeling import *
from ..rainbows import *


def test_starry_injector(signal_to_noise=10):
    residuals = SimulatedRainbow(
        signal_to_noise=signal_to_noise, tlim=[-0.6, 0.6] * u.day
    )
    injector = StarryModelInjector()
    model = injector(residuals)
    data = model * residuals
    plt.figure(figsize=(8, 3), dpi=300)
    MultiRainbow([model, residuals, data], names=["model", "noise", "data"]).imshow()
    return model, residuals, data
