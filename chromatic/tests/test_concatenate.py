from .setup_tests import *
from ..rainbows import *


def test_concatenate_simple():
    s = SimulatedRainbow()
    more_wavelengths = s.concatenate_in_wavelength(s)
    assert more_wavelengths.nwave == 2 * s.nwave
    more_times = s.concatenate_in_time(s)
    assert more_times.ntime == 2 * s.ntime
