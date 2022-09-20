from ..rainbows import *
from .setup_tests import *


def test_remove_trends():
    s = (
        SimulatedRainbow(dw=0.1 * u.micron)
        .inject_transit()
        .inject_systematics()
        .inject_noise(signal_to_noise=300)
    )
    fi, ax = plt.subplots(
        2,
        6,
        figsize=(12, 6),
        sharey="row",
        sharex=True,
        dpi=300,
        constrained_layout=True,
    )
    imkw = dict(vmin=0.98, vmax=1.02, xaxis="wavelength", colorbar=False)
    s.imshow(ax=ax[0, 0], **imkw)
    s.plot_noise_comparison(ax=ax[1, 0])
    for i, method in enumerate(
        ["median_filter", "savgol_filter", "differences", "polyfit"]
    ):
        x = s.remove_trends(method=method)
        x.imshow(ax=ax[0, 1 + i], **imkw)
        ax[0, i + 1].set_title(method)
        x.plot_noise_comparison(ax=ax[1, i + 1])
    x = s.remove_trends(method="custom", model=s.planet_model)
    ax[0, -1].set_title("custom")

    x.imshow(ax=ax[0, -1], **imkw)
    x.plot_noise_comparison(ax=ax[1, -1])
    plt.ylim(0, 0.02)
    plt.savefig(os.path.join(test_directory, "test-remove_trends.png"))
