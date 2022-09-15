from ..rainbows import *
from .setup_tests import *
from ..rainbows.helpers.history import represent_as_copypasteable


def test_history():
    np.random.seed(0)
    x = SimulatedRainbow().inject_transit().inject_noise().bin(R=5).normalize()
    h = x.history()

    for k in ["Rainbow", "inject_transit", "bin", "normalize"]:
        assert k in h

    np.random.seed(0)
    new = eval(h)
    assert x == new


def test_history_with_slicing_and_addition():
    a = SimulatedRainbow().inject_transit().normalize()[:10:2, :]
    b = a * 2
    c = a + b
    d = eval(c.history())
    assert c == d
    assert d == a * 3


def test_represent_as_copypasteable():
    for x_in in [
        np.array([1, 2]) * u.day,
        [1, 2],
        [1, 2] * u.day,
        1 * u.day,
        np.linspace(0, 1, 10),
    ]:
        string = represent_as_copypasteable(x_in)
        x_out = eval(string)
        match = np.all(np.isclose(x_in, x_out))
        print(f"{x_in} == {x_out}? {match}")
        assert match
