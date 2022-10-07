from ..writers import *


def save(self, filepath="test.rainbow.npy", format=None, **kw):
    """
    Save this `Rainbow` out to a file.

    Parameters
    ----------
    filepath : str
        The filepath pointing to the file to be written.
        (For now, it needs a `.rainbow.npy` extension.)
    format : str, optional
        The file format of the file to be written. If `None`,
        the format will be guessed automatically from the
        filepath.
    **kw : dict, optional
        All other keywords will be passed to the writer.
    """
    # figure out the best writer
    writer = guess_writer(filepath, format=format)

    # use that writer to save the file
    writer(self, filepath, **kw)
