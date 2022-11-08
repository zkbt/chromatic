from ....imports import *

__all__ = ["plot_histogram"]


def plot_histogram(
    self,
    i_wavelength,
    ax=None,
    quantity="flux",
    offset=0,
    scaling=1,
    color="auto",
    cmap=None,
    vmin=None,
    vmax=None,
    orientation="vertical",
    expected=False,
    expected_plotkw={},
    **kw,
):
    """
    Plots a histogram of the flux value for one wavelength of a rainbow.

    This makes a histogram of flux values for one wavelength, with flux bin
    values on the x-axis and the number of occurrences (= histogram) of
    flux values within a particular bin on the y-axis.

    Parameters
    ----------
    i_wavelength : integer
        The wavelength row we want to plot.
    ax : Axes, optional
        The axes into which the plot should be drawn.
    quantity : str, optional
        The quantity for which we want a histogram.
        Currently available options are ['flux', 'residuals']
    offset : float, optional
        An offset to add to each value (needed for `plot_with_model_and_residuals`)
    color : str, optional
        The color for the histogram. If 'auto', guess from wavelength.
    expected : bool, optional
        Allows user to choose whether an expected normal distribution for the
        data will be plotted over the histogram.
    expected_plotkw : dict, optional
        Keywords to pass to `plt.plot` for plotting the expected distribution.
    **kw : dictionary
        All additional keywords will be passed to `plt.hist`.
        Please see the documentation for that function for options.
        (https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html)
    """

    if ax is None:
        ax = plt.subplot()
    plt.sca(ax)

    # importing norm so we can plot normal distribution over histogram
    from scipy.stats import norm

    # pick the color of the histogram based on the wavelength
    if color == "auto":
        self._make_sure_cmap_is_defined()
        color = self.get_wavelength_color(self.wavelength[i_wavelength])

    # get the quantity for which we want to plot the histogram
    assert quantity in ["flux", "residuals", "residuals_plus_one"]
    time, flux, uncertainty = self.get_ok_data_for_wavelength(i_wavelength, y=quantity)

    # plotting histogram of row 'i' of data (wavelength 'i')
    histkw = dict(alpha=0.5)
    histkw.update(**kw)
    plt.hist(
        flux * scaling + offset,
        color=color,
        density=True,
        orientation=orientation,
        **histkw,
    )
    if orientation == "vertical":
        plt.xlabel(f"{quantity}")
        plt.ylabel(f"P({quantity})")
    elif orientation == "horizontal":
        plt.ylabel(f"{quantity}")
        plt.xlabel(f"P({quantity})")

    # option to plot expected normal distribution over histogram
    if expected == True:
        # mu (middle of the normal distribution) corresponds to the middle
        # of the flux array
        if quantity == "residuals":
            mu = 0
        else:
            mu = np.median(flux)

        # expected std corresponds to the middle of our true uncertainty
        std = np.median(uncertainty)

        # setting min/max values for scale of normal distribution
        nsigma = 4
        xmin = np.minimum(np.min(flux), mu - nsigma * std)
        xmax = np.maximum(np.max(flux), mu + nsigma * std)

        x = np.linspace(xmin, xmax, 500)
        p = norm.pdf(x, mu, std)
        plotkw = dict(color=color, alpha=0.5)
        plotkw.update(**expected_plotkw)
        if orientation == "vertical":
            plt.plot(x * scaling + offset, p / scaling, **plotkw)
        elif orientation == "horizontal":
            plt.plot(p / scaling, x * scaling + offset, **plotkw)
