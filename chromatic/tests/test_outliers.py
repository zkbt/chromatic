from ..rainbows import *
from .setup_tests import *


def test_flag_outliers():
    noiseless = SimulatedRainbow(R=10, dt=10 * u.minute).inject_transit()
    noiseless_outliers = noiseless.inject_outliers()
    noiseless_flagged = noiseless_outliers.flag_outliers()

    noisy = SimulatedRainbow(R=10, dt=10 * u.minute).inject_transit().inject_noise()
    noisy_outliers = noisy.inject_outliers()
    noisy_flagged = noisy_outliers.flag_outliers()

    kw = dict(vmin=0.9, vmax=1.1)
    fi, ax = plt.subplots(2, 3, figsize=(10, 4), dpi=300)
    noisy.imshow(ax=ax[0, 0], **kw)
    noisy_outliers.imshow(ax=ax[0, 1], **kw)
    noisy_flagged.imshow(ax=ax[0, 2], **kw)

    noiseless.imshow(ax=ax[1, 0], **kw)
    noiseless_outliers.imshow(ax=ax[1, 1], **kw)
    noiseless_flagged.imshow(ax=ax[1, 2], **kw)

    plt.savefig(os.path.join(test_directory, "demonstration-of-flagging-outliers.pdf"))
    assert np.all(
        (noisy_flagged.fluxlike["injected_outliers"] != 0)
        == noisy_flagged.fluxlike["flagged_as_outlier"]
    )
