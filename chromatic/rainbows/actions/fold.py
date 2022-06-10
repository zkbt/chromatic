from ...imports import *


def fold(self, period=None, t0=None, event="Mid-Transit"):
    """
    Fold this Rainbow to a period and reference epoch.

    Parameters
    ----------
    period : u.Quantity
        The orbital period of the planet (with astropy units of time).
    t0 : u.Quantity
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
        The folded Rainbow
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("fold", locals())

    # warn
    if (period is None) or (t0 is None):
        message = """
        Folding to a transit period requires both
        `period` and `t0` be specified. Please try again.
        """
        warnings.warn(message)
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
