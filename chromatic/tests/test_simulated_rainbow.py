from ..rainbows import *
from .setup_tests import *

def test_simulated_rainbow():

    a = SimulatedRainbow()
    b = SimulatedRainbow(tlim=[-1, 1] * u.hour, dt=1 * u.minute)
    c = SimulatedRainbow(time=np.linspace(-1, 1) * u.hour)
    d = SimulatedRainbow(wlim=[0.1, 1] * u.micron, R=100)
    assert d.wscale == "log"
    e = SimulatedRainbow(wlim=[0.1, 1] * u.micron, dw=1 * u.nm)
    assert e.wscale == "linear"
    f = SimulatedRainbow(wavelength=np.logspace(0, 1) * u.micron)
