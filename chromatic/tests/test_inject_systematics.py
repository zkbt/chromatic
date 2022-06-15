from ..rainbows import *
from .setup_tests import *


def test_inject_systematics():
    SimulatedRainbow().inject_transit().inject_noise().inject_systematics().bin(
        R=10, dt=10 * u.minute
    ).imshow_quantities()
    plt.savefig(os.path.join(test_directory, "inject_systematics_demo.pdf"))
