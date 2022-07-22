from ....imports import *

__all__ = ["plot_quantities"]


def plot_quantities(
    self,
    quantities=None,
    xaxis="time",
    maxcol=1,
    panel_size=(6, 2),
    what_is_x="index",
    filename=None,
    **kw,
):
    """
    Plot {xaxis}-like quantities as a function of {xaxis} index
    or any other {xaxis} quantity (such as "time" or "wavelength").

    Parameters
    ----------
    quantities : list like
        The {xaxis}-like quantity to plot.
    xaxis : string
        Whether the quantities are alike to 'time' or 'wave'. Default is 'time'. (Optional)
    maxcol : int
        The maximum number of columns to show (Optional).
    panel_size : tuple
        The size in inches dedicated to each panel, default is (6,2). (Optional)
    what_is_x : string
        The quantity to plot on the what_is_x, default is index. (Optional)

    """
    # decide which dictionary to plot
    if xaxis not in ["time", "wave", "wavelength"]:
        raise Exception("Unknown xaxis. Choose from [time, wave]")
    elif xaxis == "time":
        like_dict = self.timelike
    else:
        like_dict = self.wavelike

    # decide which quantities to plot
    if quantities is None:
        allkeys = like_dict.keys()
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
        constrained_layout=True,
    )
    # make the axes easier to index
    if len(allkeys) > 1:
        ax = axes.flatten()
    else:
        ax = [axes]

    # set what_is_x variable
    if what_is_x.lower() == "index":
        x_values = np.arange(0, len(like_dict[list(like_dict.keys())[0]]))
        xlab = f"{xaxis} Index"
    else:
        if what_is_x in like_dict.keys():
            x_values = like_dict[what_is_x]
            xlab = what_is_x
        else:
            raise Exception("Desired `what_is_x` quantity is not in given dictionary")

    # display each quantity
    for k, key in enumerate(allkeys):
        # make the plot (or an empty box)
        if key in like_dict.keys():
            ax[k].plot(
                x_values, like_dict[key], color=plt.cm.viridis(k / len(allkeys)), **kw
            )
            ax[k].set_xlabel(xlab)
            ax[k].set_ylabel(
                f"{key} ({u.Quantity(like_dict[key]).unit.to_string('latex_inline')})"
            )
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

    # hide any additional axes
    if k + 1 <= len(ax):
        for axi in ax[k + 1 :]:
            axi.axis("Off")
    if filename is not None:
        self.savefig(filename)
