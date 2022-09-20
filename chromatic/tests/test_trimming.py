from ..rainbows import *
from .setup_tests import *


def test_trim():
    r = SimulatedRainbow().inject_noise()
    r.fluxlike["flux"][:3, :] = np.nan
    r.fluxlike["flux"][:, -4:] = np.nan
    original_shape = r.shape
    t = r.trim()
    new_shape = t.shape
    assert new_shape[0] == original_shape[0] - 3
    assert new_shape[1] == original_shape[1] - 4
    print(r, t)
