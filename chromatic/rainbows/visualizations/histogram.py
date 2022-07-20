from ...imports import *

__all__ = ["plot_histogram"]


def plot_histogram(self, i_wavelength=0, expected=False, color="auto", **kw):
    """
    Plots a histogram of the flux value for one wavelength of a rainbow.

    Parameters
    ----------
    i_wavelength : integer
        The wavelength row we want to plot.

    expected : Boolean
        Allows user to choose whether an expected normal distribution for the
        data will be plotted over the histogram.

    kw : dictionary
        All additional keywords will be passed to `plt.hist`.
        Please see the documentation for that function for options.
        (https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html)

    Returns
    -------
        A plotted histogram of flux values for one wavelength, with flux bin
        values on the x-axis and the number of occurrences (= histogram) of
        flux values within a particular bin on the y-axis.

    """
    # importing norm so we can plot normal distribution over histogram
    from scipy.stats import norm

    # pick the color of the histogram based on the wavelength
    if color == "auto":
        self._make_sure_cmap_is_defined()
        color = self.get_wavelength_color(self.wavelength[i_wavelength])

    # get the quantity for which we want to plot the histogram
    time, flux, uncertainty = self.get_ok_data_for_wavelength(i_wavelength)

    # plotting histogram of row 'i' of data (wavelength 'i')
    plt.hist(flux, color=color, density=True, **kw)
    plt.xlabel("Flux")
    plt.ylabel("Histogram of Flux Values")

    # option to plot expected normal distribution over histogram
    if expected == True:
        # mu (middle of the normal distribution) corresponds to the middle
        # of the flux array
        mu = np.median(flux)
        # expected std corresponds to the middle of our true uncertainty
        std = np.median(uncertainty)

        # setting min/max values for scale of normal distribution
        xmin = np.min(flux)
        xmax = np.max(flux)

        x = np.linspace(xmin, xmax, 500)
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, "k", linewidth=2)
