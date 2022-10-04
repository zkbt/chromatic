from ..rainbows import *
from .setup_tests import *


def test_lightcurve_and_spectrum_summaries():
    fi, ax = plt.subplots(1, 2, sharey=True)
    s = SimulatedRainbow().inject_noise().inject_transit()
    ax[0].plot(s.time, s.get_average_lightcurve())
    ax[1].plot(s.wavelength, s.get_average_spectrum())
    plt.savefig(
        os.path.join(
            test_directory, "demonstration-of-average-lightcurve-and-spectrum.pdf"
        )
    )


def test_measured_scatter_summary():
    # make a fake Rainbow with known noise
    signal_to_noise = 100
    s = SimulatedRainbow().inject_noise(signal_to_noise=signal_to_noise)

    # mathematically, what do we expect?
    expected_sigma = 1 / signal_to_noise
    uncertainty_on_expected_sigma = expected_sigma / np.sqrt(2 * (s.ntime - 1))
    # (see equation 3.48 of Sivia and Skilling for the uncertainty on sigma)

    # loop through methods
    methods = ["standard-deviation", "MAD"]
    fi, ax = plt.subplots(
        1,
        len(methods),
        sharey=True,
        sharex=True,
        figsize=(8, 3),
        dpi=300,
        constrained_layout=True,
        facecolor="white",
    )
    for i, method in enumerate(methods):
        # calculate the measured scatter
        measured_scatter = s.get_measured_scatter(method=method)

        # plot the measured scatters + 1-sigma expectation for its uncertainty
        plt.sca(ax[i])
        plt.plot(s.wavelength, measured_scatter, marker="o", color="black")
        plt.axhline(expected_sigma, color="black", alpha=0.3)
        plt.axhspan(
            expected_sigma - uncertainty_on_expected_sigma,
            expected_sigma + uncertainty_on_expected_sigma,
            color="gray",
            alpha=0.3,
        )

        # plot range that shouldn't be exceeded in this test
        many_sigma = 10
        plt.axhspan(
            expected_sigma - uncertainty_on_expected_sigma * many_sigma,
            expected_sigma + uncertainty_on_expected_sigma * many_sigma,
            color="red",
            alpha=0.05,
            zorder=-10,
        )
        plt.ylim(0, None)
        plt.xscale("log")
        plt.xlabel(f'Wavelength ({s.wavelength.unit.to_string("latex_inline")})')
        plt.ylabel(f"Measured Scatter\n('{method}')")

        # throw an error if the measured scatter gets too large
        assert np.all(
            np.abs(measured_scatter - expected_sigma)
            < uncertainty_on_expected_sigma * many_sigma
        )

    plt.savefig(os.path.join(test_directory, "demonstration-of-measured-scatter.pdf"))
