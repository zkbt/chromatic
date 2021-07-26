from ..rainbows import *
from .setup_tests import *


def test_rainbownpy():
    filename = os.path.join(test_directory, "test.rainbow.npy")
    a = SimulatedRainbow()
    a.save(filename)
    b = Rainbow(filename)
    assert a == b
