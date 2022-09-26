from ..rainbows import *
from .setup_tests import *
import pytest


def test_rainbow_operations():
    nw, nt = 20, 40
    a = Rainbow(
        wavelength=np.linspace(0.5, 5, nw) * u.micron,
        time=np.linspace(-1, 1, nt) * u.hour,
        flux=np.ones((nw, nt)),
    )

    b = Rainbow(
        wavelength=np.linspace(0.5, 5, nw) * u.micron,
        time=np.linspace(-1, 1, nt) * u.hour,
        flux=np.zeros((nw, nt)) + 0.1,
    )

    wl_like = np.linspace(-0.01, 0.01, nw)
    t_like = np.linspace(-0.01, 0.01, nt)
    fx_like = np.ones(a.shape)

    assert (a + b).fluxlike["flux"][0][0] == 1.1
    assert (a + wl_like).fluxlike["flux"][0][0] == 0.99
    assert (a + t_like).fluxlike["flux"][0][0] == 0.99
    assert (a + fx_like).fluxlike["flux"][0][0] == 2
    assert (a + 1).fluxlike["flux"][0][0] == 2

    assert (a - b).fluxlike["flux"][0][0] == 0.9
    assert (a - wl_like).fluxlike["flux"][0][0] == 1.01
    assert (a - t_like).fluxlike["flux"][0][0] == 1.01
    assert (a - fx_like).fluxlike["flux"][0][0] == 0
    assert (a - 1).fluxlike["flux"][0][0] == 0

    assert (a * b).fluxlike["flux"][0][0] == 0.1
    assert (a * wl_like).fluxlike["flux"][0][0] == -0.01
    assert (a * t_like).fluxlike["flux"][0][0] == -0.01
    assert (a * fx_like).fluxlike["flux"][0][0] == 1
    assert (a * 1).fluxlike["flux"][0][0] == 1

    assert (a / b).fluxlike["flux"][0][0] == 1 / 0.1
    assert (a / wl_like).fluxlike["flux"][0][0] == 1 / -0.01
    assert (a / t_like).fluxlike["flux"][0][0] == 1 / -0.01
    assert (a / fx_like).fluxlike["flux"][0][0] == 1
    assert (a / 1).fluxlike["flux"][0][0] == 1

    # make sure we raise an error if it's not obvious whether we're doing wavelength or time
    with pytest.warns(match="reconsider letting them have the same size"):
        c = Rainbow(
            wavelength=np.linspace(0.5, 5, nw) * u.micron,
            time=np.linspace(-1, 1, nw) * u.hour,
            flux=np.ones((nw, nw)),
        )
        with pytest.raises(Exception):
            c * wl_like
        with pytest.raises(Exception):
            c / wl_like
        with pytest.raises(Exception):
            c + wl_like
        with pytest.raises(Exception):
            c - wl_like


def test_operations_with_uncertainty():
    a = SimulatedRainbow().inject_noise(signal_to_noise=100) + 1
    b = SimulatedRainbow().inject_noise(signal_to_noise=50) * 0.5
    for x in [
        "a",
        "b",
        "(a+0)",
        "(a+1)",
        "(a+b)",
        "(a-b)",
        "(a*b)",
        "(a/b)",
        "(a+b)/(a-b)*2",
    ]:
        print(f"{x:^38}")
        r = eval(x)

        print(f"       mean(flux) = {np.mean(r.flux)}")
        print(f"      mean(model) = {np.mean(r.model)}")
        print(f"mean(uncertainty) = {np.mean(r.uncertainty)}")
        print()

        scaled_sigma = np.std(r.residuals / r.uncertainty)
        assert np.isclose(scaled_sigma, 1, rtol=0.02)
