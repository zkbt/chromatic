from ..rainbows import *
from .setup_tests import *


def test_custom_units():
    u.Unit("(MJy/sr)^2")
    u.Unit("(DN/s)^2")
    u.Unit("microns")
    u.Unit("electrons")
