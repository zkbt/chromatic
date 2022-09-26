from ..imports import *
from ..resampling import *
from ..resampling import leftright_to_edges
from .setup_tests import *


def test_resampling():

    resample_directory = os.path.join(test_directory, "core-resample-examples")
    try:
        os.mkdir(resample_directory)
    except FileExistsError:
        pass

    N = 11
    x_edges = np.linspace(2, 4, N)
    x = (x_edges[1:] + x_edges[:-1]) / 2
    y = np.random.normal(x, 0.5)
    N_new = 29
    x_new_edges = np.linspace(1.6, 4.3, N_new)
    x_new = (x_new_edges[1:] + x_new_edges[:-1]) / 2

    a = resample_while_conserving_flux(xin=x, yin=y, xout=x_new, visualize=True)
    plt.savefig(
        os.path.join(
            resample_directory, "demonstration-of-resampling-centers-in-centers-out.pdf"
        )
    )

    b = resample_while_conserving_flux(
        yin=y, xout=x_new, xin_edges=x_edges, visualize=True
    )
    plt.savefig(
        os.path.join(
            resample_directory, "demonstration-of-resampling-edges-in-centers-out.pdf"
        )
    )

    c = resample_while_conserving_flux(
        yin=y, xin_edges=x_edges, xout_edges=x_new_edges, visualize=True
    )
    plt.savefig(
        os.path.join(
            resample_directory, "demonstration-of-resampling-edges-in-edges-out.pdf"
        )
    )

    d = resample_while_conserving_flux(
        yin=y, xin=x, xout_edges=x_new_edges, visualize=True
    )
    plt.savefig(
        os.path.join(
            resample_directory, "demonstration-of-resampling-centers-in-edges-out.pdf"
        )
    )

    for x in [b, c, d]:
        assert np.isclose(np.sum(a["y"]), np.sum(x["y"]))
        assert np.all(np.isclose(a["x_edge_lower"], x["x_edge_lower"]))
        assert np.all(np.isclose(a["x_edge_upper"], x["x_edge_upper"]))


def test_bintogrid(N_original=23):

    bin_directory = os.path.join(test_directory, "binning-examples")
    try:
        os.mkdir(bin_directory)
    except FileExistsError:
        pass

    def save_binning_example_figure():
        f = plt.gca().get_title()
        f = "demonstration-of-" + f.replace(" ", "") + ".pdf"
        f = f.replace("$\sigma$", "uncertainty")
        f = f.replace(",", "-")
        plt.savefig(os.path.join(bin_directory, f))

    for input_grid in ["uniform", "irregular"]:
        if input_grid == "uniform":
            x = np.linspace(-3, 3, N_original)
        elif input_grid == "irregular":
            x = np.sort(np.random.uniform(-3, 3, N_original))

        for uncertainties in ["no", "with"]:
            plt.close("all")
            if uncertainties == "no":
                u = None
                y = np.exp(-0.5 * x**2)
            else:
                u = np.ones(N_original) * 0.2
                y = np.random.normal(np.exp(-0.5 * x**2), u)

            for how in ["bigger", "smaller"]:
                if how == "bigger":
                    dx = np.median(np.diff(x)) * np.random.uniform(3, 5)
                elif how == "smaller":
                    dx = np.median(np.diff(x)) / np.random.uniform(3, 5)

                label = f"{input_grid} input grid, {uncertainties} $\sigma$, new bins are {how}"

                # do different options give same answer?
                a = bintogrid(x, y, unc=u, dx=dx, visualize=True)
                plt.title(f"bin to dx={dx:.3}, {label}")
                save_binning_example_figure()

                centers = a["x"]
                b = bintogrid(x, y, unc=u, newx=centers, visualize=True)
                plt.title(f"bin to {len(centers)} newx, {label}")
                save_binning_example_figure()
                assert np.all(np.isclose(a["y"], b["y"]))

                edges = leftright_to_edges(a["x_edge_lower"], a["x_edge_upper"])
                c = bintogrid(x, y, unc=u, newx_edges=edges, visualize=True)
                plt.title(f"bin to {len(edges)} newx_edges, {label}")
                save_binning_example_figure()
                assert np.all(np.isclose(a["y"], c["y"]))

            label = f"{input_grid} input grid, {uncertainties} " + r"$\sigma$"

            # does binning to self give self?
            d = bintogrid(x, y, unc=u, newx=x, visualize=True)
            plt.title(f"bin to self, {label}")
            save_binning_example_figure()
            assert np.all(np.isclose(d["x"], x))
            assert np.all(np.isclose(d["y"], y))

            # does binning by "N points" work, and match other modes?
            for nx in [1, 2, 3, 5]:
                e = bintogrid(x, y, unc=u, nx=nx, visualize=True)
                plt.title(f"bin by {nx} bins, {label}")
                save_binning_example_figure()
    plt.close("all")
