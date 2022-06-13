from ...imports import *
from ..multi import *


def compare(self, rainbows):
    """
    Compare this Rainbow to others.

    Parameters
    ----------
    rainbows : list
        A list containing one or more other Rainbow objects.
        If you only want to compare with one other Rainbow,
        supply it in a 1-element list like `.compare([other])`
    """
    try:
        rainbows.remove(self)
    except (ValueError, IndexError):
        pass
    return compare_rainbows([self] + rainbows)
