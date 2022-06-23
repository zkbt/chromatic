from ..rainbows import *
from .setup_tests import *


def test_rainbow_npy():
    filename = os.path.join(test_directory, "test.rainbow.npy")
    a = SimulatedRainbow().inject_noise()
    a.save(filename)
    b = Rainbow(filename)
    assert a == b


def test_rainbow_FITS():
    filename = os.path.join(test_directory, "test.rainbow.fits")
    a = SimulatedRainbow().inject_noise()
    a.save(filename, overwrite=True)
    b = Rainbow(filename)
    assert a == b


def test_text():
    filename = os.path.join(test_directory, "test.rainbow.txt")
    a = SimulatedRainbow().inject_noise()
    a.save(filename, overwrite=True)
    b = Rainbow(filename)
    # assert a == b


def test_guess_readers():

    assert guess_reader("some-neat-file.rainbow.npy") == from_rainbow_npy
    assert guess_reader("some-neat-file-x1dints.fits") == from_x1dints

    # eureka readers
    assert guess_reader("S3_wasp39b_ap6_bg7_SpecData.h5") == from_eureka_S3
    assert guess_reader("S4_wasp39b_ap6_bg7_LCData.h5") == from_eureka_S4
    assert (
        guess_reader(
            [
                "S5_wasp39b_ap6_bg7_Table_Save_ch0.txt",
                "S5_wasp39b_ap6_bg7_Table_Save_ch1.txt",
            ]
        )
        == from_eureka_S5
    )


def test_guess_writers():

    assert guess_writer("some-neat-file.rainbow.npy") == to_rainbow_npy
