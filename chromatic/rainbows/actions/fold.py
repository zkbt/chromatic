from ...imports import *

__all__ = ["fold", "mask_transit"]


def fold(self, period=None, t0=None, event="Mid-Transit"):
    """
    Fold this `Rainbow` to a period and reference epoch.

    This changes the times from some original time into
    a phased time, for example the time within an orbital
    period, relative to the time of mid-transit. This
    is mostly a convenience function for plotting data
    relative to mid-transit and/or trimming data based
    on orbital phase.

    Parameters
    ----------
    period : Quantity
        The orbital period of the planet (with astropy units of time).
    t0 : Quantity
        Any mid-transit epoch (with astropy units of time).
    event : str
        A description of the event that happens periodically.
        For example, you might want to switch this to
        'Mid-Eclipse' (as well as offsetting the `t0` by the
        appropriate amount relative to transit). This description
        may be used in plot labels.

    Returns
    -------
    folded : Rainbow
        The folded `Rainbow`.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("fold", locals())

    # warn
    if (period is None) or (t0 is None):
        message = """
        Folding to a transit period requires both
        `period` and `t0` be specified. Please try again.
        """
        cheerfully_suggest(message)
        return self

    # create a copy of the existing rainbow
    new = self._create_copy()

    # calculate predicted time from transit
    new.time = (((self.time - t0) + 0.5 * period) % period) - 0.5 * period
    # (the nudge by 0.5 period is to center on -period/2 to period/2)

    # change the default time label
    new.metadata["time_label"] = f"Time from {event}"

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    return new


def mask_transit(self, period, t0, duration, event="Mid-Transit"):
    """
    Mask a transit in this `Rainbow`, to focus on out-of-transit.

    Parameters
    ----------
    period : Quantity
        The orbital period of the planet (with astropy units of time).
    t0 : Quantity
        Any mid-transit epoch (with astropy units of time).
    duration : Quantity
        The total duration of the transit to remove.
    event : str
        A description of the event that happens periodically.
        For example, you might want to switch this to
        'Mid-Eclipse' (as well as offsetting the `t0` by the
        appropriate amount relative to transit). This description
        may be used in plot labels.

    Returns
    -------
    masked : Rainbow
        The `Rainbow` with the transit masked as not `ok`.
    """
    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("mask_transit", locals())

    if "fold" in self.history():
        raise ValueError(
            f"""
        It looks like this Rainbow has already been folded.
        Please run `.mask_transit` only on unfolded data.
        (It will refold while masking.)
        """
        )

    # fold and mask transit
    folded = self.fold(period=period, t0=t0, event=event)
    folded.timelike["ok"] = np.abs(folded.time) >= duration / 2

    # append the history entry to the new Rainbow
    folded._record_history_entry(h)

    return folded
