from ..rainbows import *
from .setup_tests import *


def test_ok_rows_and_columns():
    s = SimulatedRainbow().inject_noise()
    s.wavelike["ok"] = np.arange(s.nwave) > 10
    s.timelike["ok"] = np.arange(s.ntime) > 5
    s.fluxlike["ok"] = np.ones(s.shape)
    assert s.ok[0, 0] == 0
    assert s.ok[-1, -1] == 1
    assert np.all(s.ok[:10, :] == 0)
    assert np.all(s.ok[:, :5] == 0)


def test_bin_with_not_ok_data():
    for ok_fraction in [0, 0.01, 0.5, 0.99, 1]:
        a = (
            SimulatedRainbow(dt=2 * u.minute, dw=0.2 * u.micron)
            .inject_transit()
            .inject_noise()
        )
        a.ok = np.random.uniform(size=a.shape) < ok_fraction

        if ok_fraction == 0:
            with pytest.raises(RuntimeError):
                should_fal = a.bin(dw=0.7 * u.micron, dt=20 * u.minute, ok_threshold=1)
            continue

        cautious = a.bin(dw=0.7 * u.micron, dt=20 * u.minute, ok_threshold=1)
        carefree = a.bin(dw=0.7 * u.micron, dt=20 * u.minute, ok_threshold=0)
        assert np.all((cautious.ok == 1) | (cautious.ok == 0))
        if np.any(a.ok == 0):
            assert np.any((carefree.ok != 1) & (carefree.ok != 0))
