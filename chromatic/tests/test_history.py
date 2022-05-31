from ..rainbows import *
from .setup_tests import *


def test_history():
    x = SimulatedRainbow().inject_transit().bin(R=5).normalize()
    h = x.history("string")

    for k in ["Rainbow", "inject_transit", "bin", "normalize"]:
        assert k in h
