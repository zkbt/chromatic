from ..rainbows import *
from .setup_tests import *


def test_shift():
    N = 1000
    w = np.linspace(0.6, 0.7, N) * u.micron
    f = u.Quantity(np.zeros(N))
    for i in range(3):
        w0 = np.random.uniform(0.6, 0.7) * u.micron
        sigma = np.random.uniform(0.001, 0.01) * u.micron
        f += np.exp(-0.5 * ((w - w0) / sigma) ** 2)

    for v in [0, -1e4, 1e4] * u.km / u.s:
        print(f"Testing for v={v}")
        unshifted = SimulatedRainbow(wavelength=w, star_flux=f).inject_noise()
        shifted = unshifted.shift(v)
        shifted_and_then_shifted_back = shifted.shift(-v)

        fi, ax = plt.subplots(
            1, 3, sharex=True, sharey=True, figsize=(8, 4), constrained_layout=True
        )
        unshifted.imshow(ax=ax[0])
        shifted.imshow(ax=ax[1])
        shifted_and_then_shifted_back.imshow(ax=ax[2])
        plt.savefig(
            os.path.join(
                test_directory,
                "demonstration-of-shifting-wavelengths.pdf",
            )
        )

        assert np.all(
            np.isclose(unshifted.wavelength, shifted_and_then_shifted_back.wavelength)
        )

        if v.to_value("km/s") > 0:
            print(f"The velocity {v} should cause a redshift.")
            assert np.all(shifted.wavelength > unshifted.wavelength)
        elif v.to_value("km/s") == 0:
            print(f"The velocity {v} should cause no shift.")
            print("   It does!")
        elif v.to_value("km/s") < 0:
            print(f"The velocity {v} should cause a blueshift.")
            assert np.all(shifted.wavelength < unshifted.wavelength)
        print()
