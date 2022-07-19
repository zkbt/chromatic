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
    from scipy.stats import norm

    # pick the color based on the wavelength
    if color == "auto":
        self._make_sure_cmap_is_defined()
        color = self.get_wavelength_color(self.wavelength[i_wavelength])

    # get the quantity for which we want to plot the histogram
    x, y, sigma = self.get_ok_data_for_wavelength(i_wavelength)

    # option to plot normal distribution over histogram
    # This doesn't really work right now but I'm working on it.  I need to
    # straighten out 'x' because it's causing problems with scaling of the plot
    if expected == True:
        warnings.warn(f"The `expected=True` option doesn't quite work yet.")
        mu, std = norm.fit(new.flux)

        xmin, xmax = plt.xlim()
        x = np.arange(0.98, 1.08, 0.1)
        #x = np.linspace(xmin, xmax)
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, "k", linewidth=2)

    # plotting histogram of row 'i' of data (wavelength 'i')
    plt.hist(y, color=color, density=True, **kw)
    plt.xlabel("Flux")
    plt.ylabel("Histogram of Flux Values")
