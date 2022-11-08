def get(self, key, default=None):
    """
    Retrieve an attribute by its string name.
    (This is a friendlier wrapper for `getattr()`).

    `r.get('flux')` is identical to `r.flux`

    This is different from indexing directly into
    a core dictionary (for example, `r.fluxlike['flux']`),
    because it can also be used to get the results of
    properties that do calculations on the fly (for example,
    `r.residuals` in the `RainbowWithModel` class).

    Parameters
    ----------
    key : str
        The name of the attribute, property, or core dictionary item to get.
    default : any, optional
        What to return if the attribute can't be found.

    Returns
    -------
    thing : any
        The thing you were trying to get. If unavailable,
        return the `default` (which by default is `None`)
    """
    try:
        return getattr(self, key)
    except AttributeError:
        return default
