from ..rainbows import *
from .setup_tests import *


def test_astropy_times():
    s = SimulatedRainbow(time=(Time.now().jd + np.linspace(0, 1)) * u.day)
    format = "jd"
    scale = "tdb"
    original_times = s.time
    astropy_times = s.get_times_as_astropy(
        format=format, scale=scale, is_barycentric=True
    )
    s.set_times_from_astropy(astropy_times, is_barycentric=True)

    assert np.all(s.time == original_times)
    assert np.all(astropy_times == s.get_times_as_astropy())
