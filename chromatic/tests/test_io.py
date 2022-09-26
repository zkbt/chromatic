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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a.save(filename, overwrite=True)
    b = read_rainbow(filename)
    assert a == b


def test_text():
    filename = os.path.join(test_directory, "test.rainbow.txt")
    a = SimulatedRainbow().inject_noise()
    a.save(filename, overwrite=True)
    b = read_rainbow(filename)
    assert a == b


def test_xarray():

    for f in [
        "stellar-spec-planettest-modetest-codechromatic-authorzkbt.xc",
        "raw-light-curves-planettest-modetest-codechromatic-authorzkbt.xc",
        "fitted-light-curves-planettest-modetest-codechromatic-authorzkbt.xc",
    ]:
        filename = os.path.join(test_directory, f)
        print(filename)
        a = SimulatedRainbow().inject_transit().inject_systematics().inject_noise()

        with pytest.warns(match="required metadata keyword"):
            a.save(filename)
        b = read_rainbow(filename)
        assert a == b


def test_rainbow_npy():
    filename = os.path.join(test_directory, "test.rainbow.npy")
    a = SimulatedRainbow().inject_noise()
    a.save(filename)
    b = Rainbow(filename)
    assert a == b


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

    # xarray common format
    assert guess_reader("stellar-spec-wow.xc") == from_xarray_stellar_spectra
    assert guess_reader("raw-light-curves-wow.xc") == from_xarray_raw_light_curves
    assert guess_reader("fitted-light-curve-wow.xc") == from_xarray_fitted_light_curves


def test_guess_writers():

    assert guess_writer("some-neat-file.rainbow.npy") == to_rainbow_npy

    # xarray common format
    assert guess_writer("stellar-spec-wow.xc") == to_xarray_stellar_spectra
    assert guess_writer("raw-light-curves-wow.xc") == to_xarray_raw_light_curves
    assert guess_writer("fitted-light-curve-wow.xc") == to_xarray_fitted_light_curves
