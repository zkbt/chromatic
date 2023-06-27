from ...imports import *

# FIXME:
# - include more careful handling of different keys in different rainbows
# - write more flexible `concatenate_flexibly` (?) to allow different wavelengths *and* times
# - write wrapper `concatenate` to choose the approproate concatenator
# - think about what to do with merged histories?

__all__ = ["concatenate_in_time", "concatenate_in_wavelength"]


def fractional_difference(a, b):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.abs((a - b) / (a + b) * 2)


def concatenate_in_time(self, other, maximum_fractional_difference=0.01):
    """
    Merge another Rainbow into this one,
    assuming their wavelengths are *exactly* the same.

    Parameters
    ----------
    other : Rainbow
       The other Rainbow to merge into this one.

    maximum_fractional_difference : float
        If the fractional difference between arrays that
        should be the same (such as wavelength) exceeds
        this amount, the user will be warned about it.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("concatenate_in_time", locals())

    # make sure wavelength grid shapes match, at least
    assert self.nwave == other.nwave

    # create an exact copy of the first object
    new = self._create_copy()

    # loop through wavelengths
    for k in self.wavelike:
        f = fractional_difference(self.wavelike[k], other.wavelike[k])
        if np.any(f > maximum_fractional_difference):
            cheerfully_suggest(
                f"""
            wavelike['{k}'] differs fractionally by a maximum of {np.max(f):.3g}
            between `self` and `other` ({np.sum(f > maximum_fractional_difference)} differences greater than {maximum_fractional_difference:.3g}).
            `self` values will be used; `other` values will be ignored.
            """
            )

    # loop through timelike quantities
    for k in self.timelike:
        new.timelike[k] = np.hstack([self.timelike[k], other.timelike[k]])

    # loop through fluxlike quantities
    for k in self.fluxlike:
        new.fluxlike[k] = np.hstack([self.fluxlike[k], other.fluxlike[k]])

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new


def concatenate_in_wavelength(self, other, maximum_fractional_difference=0.01):
    """
    Merge another Rainbow into this one,
    assuming their times are *exactly* the same.

    Parameters
    ----------
    other : Rainbow
       The other Rainbow to merge into this one.

    maximum_fractional_difference : float
        If the fractional difference between arrays that
        should be the same (such as time) exceeds
        this amount, the user will be warned about it.
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("concatenate_in_wavelength", locals())

    # make sure wavelength grid shapes match, at least
    assert self.ntime == other.ntime

    # create an exact copy of the first object
    new = self._create_copy()

    # loop through times
    for k in self.timelike:
        f = fractional_difference(self.timelike[k], other.timelike[k])
        if np.any(f > maximum_fractional_difference):
            cheerfully_suggest(
                f"""
            timelike['{k}'] differs fractionally by a maximum of {np.max(f):.3g}
            between `self` and `other` ({np.sum(f > maximum_fractional_difference)} differences greater than {maximum_fractional_difference:.3g}).
            `self` values will be used; `other` values will be ignored.
            """
            )

    # loop through wavelike quantities
    for k in self.wavelike:
        new.wavelike[k] = np.hstack([self.wavelike[k], other.wavelike[k]])

    # loop through fluxlike quantities
    for k in self.fluxlike:
        new.fluxlike[k] = np.hstack([self.fluxlike[k], other.fluxlike[k]])

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new
