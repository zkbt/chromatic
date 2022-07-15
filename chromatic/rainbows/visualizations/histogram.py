from ...imports import *

__all__ = ["plot_histogram"]


def plot_histogram(self, i=0, bins=None, rwidth=None, expected=False):
    """
    Plots a histogram of the flux value for one wavelength of a rainbow.

    Parameters
    ----------
    i : integer
        The wavelength row we want to plot.

    bins : integer
        Number of bins desired to use in the plot.  Uses default value of 10
        unless otherwise chosen by user.

    rwidth : float
        Width of the bins.  Automatically chosen unless
        user wants to change it.

    expected : Boolean
        Allows user to choose whether an expected normal distribution for the
        data will be plotted over the histogram.

    Returns
    -------
        A plotted histogram of flux values for one wavelength, with flux bin
        values on the x-axis and the number of occurrences (= histogram) of
        flux values within a particular bin on the y-axis.

    """
    from scipy.stats import norm
    # Creating copy of rainbow to use to plot histogram
    new = self._create_copy()

    # option to plot normal distribution over histogram
    # This doesn't really work right now but I'm working on it.  I need to
    # straighten out 'x' because it's causing problems with scaling of the plot
    if expected == True:
        mu, std = norm.fit(new.flux)

        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax)
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, 'k', linewidth=2)



    # plotting histogram of row 'i' of data (wavelength 'i')
    plt.hist(new.flux[i,:], bins, rwidth, color = 'red', density=True)
    plt.xlabel('Flux')
    plt.ylabel('Occurrence of Flux Values')
    plt.title('Histogram of Flux for One Wavelength')
