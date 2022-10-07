from ...imports import *

__all__ = [
    "_add_panel_labels",
    "get_animation_writer_and_displayer",
    "_scatter_timelike_or_wavelike",
    "_get_unit_string",
    "_get_plot_directory",
    "_label_plot_file",
    "savefig",
]


def get_animation_writer_and_displayer(filename="animation.html", **kw):
    """
    Create the right animation writer based on filename.

    Parameters
    ----------
    filename : str
        The filename of the movie to create.

    Returns
    -------
    writer : MovieWriter
        The matplotlib writer object.
    displayer : ?
        The
    """

    # define the options
    writers = {"html": ani.HTMLWriter, "mp4": ani.FFMpegWriter, "gif": ani.PillowWriter}
    warnings = {
        "html": "Please try `pip insall matplotlib --upgrade` and rerunning?",
        "mp4": "Please try `conda install ffmpeg` and rerunning?",
        "gif": "Please try `pip insall matplotlib --upgrade` and rerunning?",
    }
    from IPython.display import HTML, Video, Image

    displayers = {"html": HTML, "mp4": Video, "gif": Image}

    # get the writer object
    suffix = filename.split(".")[-1]
    writer = writers[suffix](**kw)
    displayer = displayers[suffix]

    if writer.isAvailable():
        return writer, displayer
    else:
        raise ValueError(
            f"""
        The writer {writer} needed for your `.{suffix}` file is not available.
        {warnings[suffix]}
        """
        )


def _add_panel_labels(axes, preset="inside", **kw):
    """
    Add (a), (b), (c) labels to a group of axes.

    Parameters
    ----------
    ax : list, AxesSubplot
        The axes into which the labels should be drawn.
    preset : str
        A few presets for where to put the labels relative to
        upper left corner of each panel. Options are ['inside', 'above']
    kw : dict
        All addition keywords will be passed to `plt.text`,
        and they will overwrite defaults.
    """

    textkw = dict(x=0.02, y=0.98, va="top", ha="left")
    if preset == "inside":
        textkw.update(x=0.02, y=0.98, va="top")
    elif preset == "outside":
        textkw.update(x=0, y=1.02, va="bottom")
    textkw.update(**kw)

    letters = "abcdefghijklmnopqrstuvwxyz"
    for i, a in enumerate(axes.flatten()):
        textkw["s"] = f"({letters[i]})"
        textkw["transform"] = a.transAxes
        a.text(**textkw)


def _scatter_timelike_or_wavelike(
    self,
    x,
    y,
    ylabel="",
    ax=None,
    color="auto",
    cmap=None,
    vmin=None,
    vmax=None,
    t_unit="day",
    w_unit="micron",
    wavelength_for_color=None,
    percentiles=(0.1, 99.9),
    ylim=(None, None),
    scatterkw={},
    **kw,
):
    """
    Scatter a quantity with either wavelength or time on the xaxis.
    (This is a helper to save code in other places.)

    Parameters
    ----------
    x : Quantity
        Either time or wavelength, with units (which will
        used to set the xlabel and/or point colors).
    y : array, Quantity
        The values to plot on the y axis.
    ylabel : string
        The ylabel for the plot. (xlabel will be guessed)
    ax : Axes
        The axes into which this plot should go.
    cmap : str, Colormap
        The color map to use for expressing wavelength.
    vmin : Quantity
        The minimum value to use for the wavelength colormap.
    vmax : Quantity
        The maximum value to use for the wavelength colormap.
    w_unit : str, Unit
        The unit for plotting wavelengths (if needed).
    t_unit : str, Unit
        The unit for plotting times (if needed).
    wavelength_for_color : Quantity
        It you're plotting a timelike quantity, and you want
        to set the color automatically based on wavelength,
        supply the wavelength for that color here.
    **kw : dict
        All additional keywords will be passed to `plt.scatter`.
    """
    with quantity_support():

        # make sure ax is set up
        if ax is None:
            ax = plt.subplot()
        plt.sca(ax)

        if x.unit.is_equivalent("m"):
            w_unit = u.Unit(w_unit)
            xlabel = f'{self._wave_label} ({w_unit.to_string("latex_inline")})'
            c = self.wavelength
        elif x.unit.is_equivalent("s"):
            t_unit = u.Unit(t_unit)
            xlabel = f'{self._time_label} ({t_unit.to_string("latex_inline")})'
            c = wavelength_for_color
        else:
            cheerfully_suggest(
                f"""
            Your requested xaxis='{xaxis} is not allowed.
            Please choose 'time' or 'wavelength'.
            """
            )
        if color == "auto":
            # make sure that the wavelength-based colormap is defined
            self._make_sure_cmap_is_defined(cmap=cmap, vmin=vmin, vmax=vmax)
            scatterkw = dict(c=c, cmap=self.cmap, norm=self.norm)
        else:
            scatterkw = dict(color=color)
        scatterkw.update(**kw)
        assert np.shape(x) == np.shape(y)

        plt.scatter(x, y, **scatterkw)
        ylim = list(ylim)
        for i in [0, 1]:
            if ylim[i] is None:
                ylim[i] = np.percentile(y, percentiles[i])
        plt.ylim(*ylim)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(self.get("title"))


def _get_unit_string(y):
    y_unit = u.Quantity(y).unit
    if y_unit == u.Unit(""):
        unit_string = "unitless"
    else:
        unit_string = y_unit.to_string("latex_inline")
    return unit_string


def _get_plot_directory(self):
    try:
        return self._plot_directory
    except AttributeError:
        directory = ""
        for option in ["title", "label", "name", "directory"]:
            x = self.get(option)
            if x is not None:
                directory = x.replace("|", "-").replace(" ", "")
                break
        self._plot_directory = directory
        if self._plot_directory != "":
            try:
                os.mkdir(self._plot_directory)
            except FileExistsError:
                pass
        return self._plot_directory


def _label_plot_file(self, filename):
    directory = self._get_plot_directory()
    if directory == "":
        return filename
    else:
        return os.path.join(directory, filename.replace(".", f"-{directory}."))


def savefig(self, filename="test.png", dpi=300, **kw):
    plt.savefig(self._label_plot_file(filename), dpi=dpi, **kw)
