"""
Define a writer for chromatic .rainbow.npy files.
"""

# import the general list of packages
from ...imports import *
from ...version import version

# define list of the only things that will show up in imports
__all__ = ["to_rainbow_npy"]


def to_rainbow_npy(self, filepath, **kw):
    """
    Write a Rainbow to a file in the .rainbow.npy format.

    Parameters
    ----------

    self : Rainbow
        The object to be saved.

    filepath : str
        The path to the file to write.
    """

    assert ".rainbow.npy" in filepath

    # populate a dictionary containing the four core dictionaries
    dictionary_to_save = self._get_core_dictionaries()

    # save that to a file
    np.save(filepath, [dictionary_to_save, version()], allow_pickle=True)
