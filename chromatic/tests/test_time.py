from ..rainbows import *
from .setup_tests import *


def test_astropy_times():
    s = SimulatedRainbow()
    format = "jd"
    scale = "tdb"
    original_times = s.time
    astropy_times = s.get_times_as_astropy(format=format, scale=scale)
    s.set_times_from_astropy(astropy_times)

    assert np.all(s.time == original_times)
    assert np.all(astropy_times == s.get_times_as_astropy())
