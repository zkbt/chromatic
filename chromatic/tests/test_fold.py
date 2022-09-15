from ..rainbows import *
from .setup_tests import *


def test_fold():
    original_time = np.linspace(-0.1, 0.1) * u.day

    for i in range(5):
        period = np.random.uniform(1, 20) * u.day
        t0 = np.random.uniform(-10, 10) * u.day
        N = np.random.randint(-10, 10)
        s = SimulatedRainbow(time=original_time + t0 + period * N).inject_noise()
        assert np.all(
            np.isclose(s.fold(period=period, t0=t0).time, original_time, atol=1e-12)
        )
