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

    This will default to None if the attribute can't be found.
    """
    try:
        return getattr(self, key)
    except AttributeError:
        return default
