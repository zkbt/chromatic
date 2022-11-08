from ..rainbows import *
from .setup_tests import *


def test_fold(N=5):

    for i in range(5):
        original_time = (
            np.linspace(np.random.uniform(-0.5, -0.2), np.random.uniform(0.2, 0.5), 500)
            * u.day
        )
        period = np.random.uniform(1, 10) * u.day
        t0 = np.random.uniform(-10, 10) * u.day
        N = np.random.randint(-10, 10)
        s = (
            SimulatedRainbow(time=original_time + t0 + period * N, R=5)
            .inject_transit(P=period.to_value("day"), t0=t0.to_value("day"))
            .inject_noise(signal_to_noise=1000)
        )
        f = s.fold(period=period, t0=t0)
        assert np.all(np.isclose(f.time, original_time, atol=1e-12))
    fi, ax = plt.subplots(1, 2, figsize=(8, 3), constrained_layout=True)
    s.imshow(ax=ax[0])
    f.imshow(ax=ax[1])
    plt.savefig(
        os.path.join(test_directory, "demonstration-of-folding-to-period-and-t0.pdf")
    )
