# is this necessary or pointless?
__all__ = ["get_flattened", "unflatten"]


def get_flattened(self, quantity="flux"):
    """
    Convert a 2D flux-like quantity into a 1D array.

    Flux-like quantities have a shape of (nwave, ntime).
    This function flattens them into a 1D array with
    a length of (nwave*time). It can be undone with
    `.unflatten()`.

    Parameters
    ----------
    quantity : str
        The name of the quantity to flatten.

    Returns
    -------
    flat : Quantity
        The flattened 1D array with size (nwave*ntime).
    """
    return self.get(quantity).flatten()


def unflatten(self, flat_array):
    """
    Convert a 1D array into a 2D flux-like quantity.

    Undo the effects of `.flatten()` by converting
    a 1D array of length (nwave*ntime) into a 2D array
    with a shape of (nwave, ntime).

    Parameters
    ----------
    flat_array : Quantity
        The flattened array.

    Returns
    -------
    unflattened : Quantity
        The 2D array, with shape (nwave, ntime).
    """
    return flat_array.reshape(self.shape)
