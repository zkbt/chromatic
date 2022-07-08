from ...imports import *

__all__ = ["_add_panel_labels", "_get_animation_writer_and_displayer"]


def _get_animation_writer_and_displayer(filename="animation.html", **kw):
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
        {warnings[k]}
        """
        )


def _add_panel_labels(axes, preset="inside", **kw):
    """
    Add (a), (b), (c) labels to a group of axes.

    Parameters
    ----------
    ax : list or array of matplotlib.axes._subplots.AxesSubplot objects
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
        textkw.update(x=0.02, y=0.98, va="top", color="white")
    elif preset == "outside":
        textkw.update(x=0, y=1.02, va="bottom", color="black")
    textkw.update(**kw)

    letters = "abcdefghijklmnopqrstuvwxyz"
    for i, a in enumerate(axes.flatten()):
        textkw["s"] = f"({letters[i]})"
        textkw["transform"] = a.transAxes
        a.text(**textkw)
