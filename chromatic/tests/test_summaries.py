from ..rainbows import *
from .setup_tests import *


def test_lightcurve_and_spectrum_summaries():
    fi, ax = plt.subplots(1, 2, sharey=True)
    s = SimulatedRainbow().inject_noise().inject_transit()
    ax[0].plot(s.time, s.get_lightcurve())
    ax[1].plot(s.wavelength, s.get_spectrum())
