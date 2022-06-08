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
