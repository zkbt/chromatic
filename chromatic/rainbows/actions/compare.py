from ...imports import *
from ..multi import *

__all__ = ["compare"]


def compare(self, rainbows):
    """
    Compare this `Rainbow` to others.

    (still in development) This connects the current `Rainbow`
    to a collection of other `Rainbow` objects, which can then
    be visualized side-by-side in a uniform way.

    Parameters
    ----------
    rainbows : list
        A list containing one or more other `Rainbow` objects.
        If you only want to compare with one other `Rainbow`,
        supply it in a 1-element list like `.compare([other])`

    Returns
    -------
    rainbow : MultiRainbow
        A `MultiRainbow` comparison object including all input `Rainbow`s
    """
    try:
        rainbows.remove(self)
    except (ValueError, IndexError):
        pass
    return compare_rainbows([self] + rainbows)
