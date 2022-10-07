from ...imports import *

__all__ = ["trim", "trim_times", "trim_wavelengths"]


def trim_times(self, just_edges=True, when_to_give_up=1, minimum_acceptable_ok=1):
    """
    Trim times that are all (or mostly) useless.

    Parameters
    ----------
    just_edges : bool, optional
        Should we only trim the outermost bad time bins?
            `True` = Just trim off the bad edges and keep
            interior bad values. Keeping interior data, even if
            they're all bad, often helps to make for more
            intuititive imshow plots.
            `False` = Trim off every bad time, whether it's on
            the edge or somewhere in the middle of the dataset.
            The resulting Rainbow will be smaller, but it might
            be a little tricky to visualize with imshow.
    when_to_give_up : float, optional
        The fraction of wavelengths that must be nan or not OK
        for the entire time to be considered bad (default = 1).
            `1.0` = trim only if all wavelengths are bad
            `0.5` = trim if more than 50% of wavelengths are bad
            `0.0` = trim if any wavelengths are bad
    minimum_acceptable_ok : float, optional
        The numbers in the `.ok` attribute express "how OK?" each
        data point is, ranging from 0 (not OK) to 1 (super OK).
        In most cases, `.ok` will be binary, but there may be times
        where it's intermediate (for example, if a bin was created
        from some data that were not OK and some that were).
        The `minimum_acceptable_ok` parameter allows you to specify what
        level of OK-ness for a point to not get trimmed.

    Returns
    -------
    trimmed : Rainbow
        The trimmed `Rainbow`.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("trim_times", locals())

    # figure out which times should be considered bad
    is_nan = np.isnan(self.flux)
    isnt_ok = self.ok < minimum_acceptable_ok
    fraction_bad = np.sum(is_nan | isnt_ok, axis=self.waveaxis) / self.nwave
    should_be_kept = fraction_bad < when_to_give_up

    # only make cuts on the edges (if desired)
    if just_edges:
        isnt_before_first = np.cumsum(should_be_kept) > 0
        isnt_after_last = (np.cumsum(should_be_kept[::-1]) > 0)[::-1]
        isnt_edge = isnt_before_first & isnt_after_last
        should_be_kept = should_be_kept | isnt_edge

    # actually try the Rainbow
    new = self[:, should_be_kept]
    new._remove_last_history_entry()

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new


def trim_wavelengths(self, just_edges=True, when_to_give_up=1, minimum_acceptable_ok=1):
    """
    Trim wavelengths that are all (or mostly) useless.

    Parameters
    ----------
    just_edges : bool, optional
        Should we only trim the outermost bad wavelength bins?
            `True` = Just trim off the bad edges and keep
            interior bad values. Keeping interior data, even if
            they're all bad, often helps to make for more
            intuititive imshow plots.
            `False` = Trim off every bad wavelength, whether it's on
            the edge or somewhere in the middle of the dataset.
            The resulting Rainbow will be smaller, but it might
            be a little tricky to visualize with imshow.
    when_to_give_up : float, optional
        The fraction of times that must be nan or not OK
        for the entire wavelength to be considered bad (default = 1).
            `1.0` = trim only if all times are bad
            `0.5` = trim if more than 50% of times are bad
            `0.0` = trim if any times are bad
    minimum_acceptable_ok : float, optional
        The numbers in the `.ok` attribute express "how OK?" each
        data point is, ranging from 0 (not OK) to 1 (super OK).
        In most cases, `.ok` will be binary, but there may be times
        where it's intermediate (for example, if a bin was created
        from some data that were not OK and some that were).
        The `minimum_acceptable_ok` parameter allows you to specify what
        level of OK-ness for a point to not get trimmed.

    Returns
    -------
    trimmed : Rainbow
        The trimmed `Rainbow`.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("trim_wavelengths", locals())

    # figure out which wavelengths should be considered bad
    is_nan = np.isnan(self.flux)
    isnt_ok = self.ok < minimum_acceptable_ok
    fraction_bad = np.sum(is_nan | isnt_ok, axis=self.timeaxis) / self.ntime
    should_be_kept = fraction_bad < when_to_give_up

    # only make cuts on the edges (if desired)
    if just_edges:
        isnt_before_first = np.cumsum(should_be_kept) > 0
        isnt_after_last = (np.cumsum(should_be_kept[::-1]) > 0)[::-1]
        isnt_edge = isnt_before_first & isnt_after_last
        should_be_kept = should_be_kept | isnt_edge

    # actually try the Rainbow
    new = self[should_be_kept, :]
    new._remove_last_history_entry()

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new


def trim(self, just_edges=True, when_to_give_up=1, minimum_acceptable_ok=1):
    """
    Trim away bad wavelengths and/or times.

    If entire wavelengths or times are marked as not `ok`,
    we can probably remove them to simplify calculations
    and visualizations. This function will trim those away,
    by default only removing problem rows/columns on the ends,
    to maintain a contiguous block.

    Parameters
    ----------
    just_edges : bool, optional
        Should we only trim the outermost bad wavelength bins?
            `True` = Just trim off the bad edges and keep
            interior bad values. Keeping interior data, even if
            they're all bad, often helps to make for more
            intuititive imshow plots.
            `False` = Trim off every bad wavelength, whether it's on
            the edge or somewhere in the middle of the dataset.
            The resulting Rainbow will be smaller, but it might
            be a little tricky to visualize with imshow.
    when_to_give_up : float, optional
        The fraction of times that must be nan or not OK
        for the entire wavelength to be considered bad (default = 1).
            `1.0` = trim only if all times are bad
            `0.5` = trim if more than 50% of times are bad
            `0.0` = trim if any times are bad
    minimum_acceptable_ok : float, optional
        The numbers in the `.ok` attribute express "how OK?" each
        data point is, ranging from 0 (not OK) to 1 (super OK).
        In most cases, `.ok` will be binary, but there may be times
        where it's intermediate (for example, if a bin was created
        from some data that were not OK and some that were).
        The `minimum_acceptable_ok` parameter allows you to specify what
        level of OK-ness for a point to not get trimmed.

    Returns
    -------
    trimmed : Rainbow
        The trimmed `Rainbow`.
    """

    trimmed = self.trim_times(
        when_to_give_up=when_to_give_up,
        just_edges=just_edges,
        minimum_acceptable_ok=minimum_acceptable_ok,
    )
    trimmed = trimmed.trim_wavelengths(
        when_to_give_up=when_to_give_up,
        just_edges=just_edges,
        minimum_acceptable_ok=minimum_acceptable_ok,
    )

    return trimmed
