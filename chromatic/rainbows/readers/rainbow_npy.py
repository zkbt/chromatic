"""
Define a reader for chromatic .rainbow.npy files.
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_rainbow_npy"]


def from_rainbow_npy(rainbow, filepath):
    """
    Populate a Rainbow from a file in the .rainbow.npy format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.

    filepath : str
        The path to the file to load, which should probably
        have an extension of `.rainbow.npy`
    """

    # read in your file, however you like
    loaded_core_dictionaries, version_used = np.load(filepath, allow_pickle=True)

    for k in rainbow._core_dictionaries:
        vars(rainbow)[k] = loaded_core_dictionaries[k]
