from ..imports import *
from ..resampling import *
from .setup_tests import *


def test_resampling():
    N = 11
    x_edges = np.linspace(3, 8, N)
    x = (x_edges[1:] + x_edges[:-1]) / 2
    y = np.random.normal(x ** 2, 0.5)
    N_new = 29
    x_new_edges = np.linspace(1, 12, N_new)
    x_new = (x_new_edges[1:] + x_new_edges[:-1]) / 2

    a = resample_while_conserving_flux(xin=x, yin=y, xout=x_new, visualize=True)
    plt.savefig(os.path.join(test_directory, "resampling-demonstration-a.pdf"))

    b = resample_while_conserving_flux(
        yin=y, xout=x_new, xin_edges=x_edges, visualize=True
    )
    plt.savefig(os.path.join(test_directory, "resampling-demonstration-b.pdf"))

    c = resample_while_conserving_flux(
        yin=y, xin_edges=x_edges, xout_edges=x_new_edges, visualize=True
    )
    plt.savefig(os.path.join(test_directory, "resampling-demonstration-c.pdf"))

    d = resample_while_conserving_flux(
        yin=y, xin=x, xout_edges=x_new_edges, visualize=True
    )
    plt.savefig(os.path.join(test_directory, "resampling-demonstration-d.pdf"))

    for x in [b, c, d]:
        assert np.sum(a) == np.sum(x)
