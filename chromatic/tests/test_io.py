from ..rainbows import *
from .setup_tests import *


def test_rainbownpy():
    filename = os.path.join(test_directory, "test.rainbow.npy")
    a = SimulatedRainbow()
    a.save(filename)
    b = Rainbow(filename)
    assert a == b


def test_guess_readers():

    assert guess_reader("some-neat-file.rainbow.npy") == from_rainbownpy
    assert guess_reader("some-neat-file-x1dints.fits") == from_x1dints
    assert guess_reader("S3_neat-planet_Save.txt") == from_eureka


def test_guess_writers():

    assert guess_writer("some-neat-file.rainbow.npy") == to_rainbownpy
