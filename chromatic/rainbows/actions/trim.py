from ...imports import *


def trim_nan_times(self, threshold=1.0):
    """
    Trim times that are all (or mostly) not numbers.

    Parameters
    ----------
    threshold
        The fraction of wavelengths that must be nan in order for
        the entire time to be considered bad (default = 1).
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("trim_nan_times", locals())

    # figure out which times are good enough to keep
    fraction_of_nans = np.sum(np.isnan(self.flux), axis=self.waveaxis) / self.nwave
    indices_of_times_to_keep = fraction_of_nans < threshold

    new = self[:, indices_of_times_to_keep]

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new


def trim_nan_wavelengths(self, threshold=1.0):
    """
    Trim wavelengths that are all (or mostly) not numbers.

    Parameters
    ----------
    threshold
        The fraction of times that must be nan in order for
        the entire wavelength to be considered bad (default = 1).
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("trim_nan_wavelengths", locals())

    # figure out which times are good enough to keep
    fraction_of_nans = np.sum(np.isnan(self.flux), axis=self.timeaxis) / self.ntime
    indices_of_wavelengths_to_keep = fraction_of_nans < threshold

    new = self[indices_of_wavelengths_to_keep, :]

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the new Rainbow
    return new


def trim(self, threshold=1.0):
    """
    Trim wavelengths or times that are all (or mostly) not numbers.

    Parameters
    ----------
    threshold
        The fraction of a particular time/wavelengths that must be nan
        in order for the entire wavelength/time to be considered bad
        (default = 1).
    """

    trimmed = self.trim_nan_times(threshold=threshold)
    trimmed = trimmed.trim_nan_wavelengths(threshold=threshold)

    return trimmed
