from ....imports import *

__all__ = ["imshow_quantities"]


def imshow_quantities(
    self, quantities=None, maxcol=3, panel_size=(5, 4), filename=None, **kw
):
    """
    imshow fluxlikes as a function of time (x = time, y = wavelength, color = flux).

    Parameters
    ----------
    quantities : None, list, optional
        The fluxlike quantity to imshow.
    maxcol : int, optional
        The maximum number of columns to show (Optional).
    panel_size : tuple, optional
        The (approximate) size of a single panel, which will
        be used to set the overall figsize based on the
        number of rows and columns (Optional).
    **kw : dict, optional
        Additional keywords will be passed on to `imshow`

    """

    # decide which quantities to plot
    if quantities is None:
        allkeys = self.fluxlike.keys()
    else:
        allkeys = quantities[:]

    # set up the geometry of the grid
    if len(allkeys) > maxcol:
        rows = int(np.ceil(len(allkeys) / maxcol))
        cols = maxcol
    else:
        rows = 1
        cols = np.min([len(allkeys), maxcol])

    # create the figure and grid of axes
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(cols * panel_size[0], rows * panel_size[1]),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )

    # make the axes easier to index
    if len(allkeys) > 1:
        ax = axes.flatten()
    else:
        ax = [axes]

    # display each quantity
    for k, key in enumerate(allkeys):
        # make the imshow (or an empty box)
        if key in self.fluxlike.keys():
            self.imshow(quantity=key, ax=ax[k], **kw)
        else:
            ax[k].text(
                0.5,
                0.5,
                f"No {key}",
                transform=ax[k].transAxes,
                ha="center",
                va="center",
            )
        # add a title for each box
        ax[k].set_title(key)

        # hide xlabel except on the bottom row
        if k < (len(allkeys) - cols):
            ax[k].set_xlabel("")
        else:
            ax[k].tick_params(labelbottom=True)

        # hide ylabel except on the left column
        if (k % cols) > 0:
            ax[k].set_ylabel("")

    # hide any additional axes
    if k + 1 <= len(ax):
        for axi in ax[k + 1 :]:
            axi.axis("Off")
    if filename is not None:
        self.savefig(filename)
