from ..rainbows import *
from .setup_tests import *


def test_rainbow_basics():

    nw, nt = 24, 48
    r = Rainbow(
        wavelength=np.linspace(1, 2, nw) * u.micron,
        time=np.linspace(-1, 1, nt) * u.hour,
        flux=np.ones((nw, nt)),
    )

    # test the basic sizes and shapes
    assert r.shape == (nw, nt)
    assert r.nwave == nw
    assert r.ntime == nt
    assert r.nflux == nw * nt

    # testing slicing, indexing, masking
    assert r[:, :].shape == (nw, nt)
    assert r[::2, ::3].shape == (nw / 2, nt / 3)
    assert r[2:4, :].shape == (2, nt)
    assert r[np.arange(2, 4), :].shape == (2, nt)
    assert r[r.wavelength < 1.5 * u.micron, :].shape == (nw / 2, nt)
    assert r[:, r.time < 0 * u.day].shape == (nw, nt / 2)

    assert r[r.wavelength < 1.5 * u.micron, r.time < 0 * u.day].shape == (
        nw / 2,
        nt / 2,
    )
