from ..imports import *


def __add__(self, object):
    """
    Add the flux of a rainbow and an input array (or another rainbow)
    and output in a new rainbow object.
    Currently only supports flux addition.

    Parameters
    ----------
    object : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow object with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow() with the same parameters as self but with added
    flux.
    """

    # Create new Rainbow() to store results in.
    result = self._create_copy()

    try:  # Object is rainbow
        if np.array_equal(
            self.wavelike["wavelength"], object.wavelike["wavelength"]
        ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
            result.fluxlike["flux"] += object.fluxlike["flux"]
        else:
            print("Objects do not share wavelength/time axes")
            return

    except AttributeError:  # Object is not rainbow

        if type(object) == float or type(object) == int:  # Float or integer.
            result.fluxlike["flux"] += object

        elif self.shape == object.shape:  # fluxlike
            result.fluxlike["flux"] += object

        elif len(object) == len(
            self.wavelike["wavelength"]
        ):  # Add along wavelength axis

            result.fluxlike["flux"] += np.transpose([object] * self.shape[1])
            if self.nwave == self.ntime:
                raise RuntimeError(
                    f"{self} has same number of wavelengths and times; we can't tell which is which."
                )

        elif len(object) == len(self.timelike["time"]):  # Add along time axis.
            result.fluxlike["flux"] += np.tile(object, (self.shape[0], 1))

        else:
            print("Invalid shape in addition: " + str(object.shape))
            return
    return result


def __sub__(self, object):
    """
    Subtract the flux of a rainbow from an input array (or another rainbow)
    and output in a new rainbow object.
    Currently only supports flux subtraction.

    Parameters
    ----------
    object : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow object with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow() with the same parameters as self but with subtracted
    flux.
    """

    # Create new Rainbow() to store results in.
    result = self._create_copy()

    try:  # Object is rainbow.
        if np.array_equal(
            self.wavelike["wavelength"], object.wavelike["wavelength"]
        ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
            result.fluxlike["flux"] -= object.fluxlike["flux"]
        else:
            print("Objects do not share wavelength/time axes")
            return

    except AttributeError:

        if type(object) == float or type(object) == int:  # Float or integer.
            result.fluxlike["flux"] -= object

        elif self.shape == object.shape:  # fluxlike
            result.fluxlike["flux"] -= object

        elif len(object) == len(
            self.wavelike["wavelength"]
        ):  # Add along wavelength axis
            result.fluxlike["flux"] -= np.transpose([object] * self.shape[1])
            if self.nwave == self.ntime:
                raise RuntimeError(
                    f"{self} has same number of wavelengths and times; we can't tell which is which."
                )

        elif len(object) == len(self.timelike["time"]):  # Add along time axis.
            result.fluxlike["flux"] -= np.tile(object, (self.shape[0], 1))

        else:
            print("Invalid shape in subtraction: " + str(object.shape))
            return
    return result


def __mul__(self, object):
    """
    Multiply the flux of a rainbow and an input array (or another rainbow)
    and output in a new rainbow object.
    Currently only supports flux multiplication.

    Parameters
    ----------
    object : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow object with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow() with the same parameters as self but with multiplied
    flux.
    """

    # Create new Rainbow() to store results in.
    result = self._create_copy()

    try:  # Object is rainbow
        if np.array_equal(
            self.wavelike["wavelength"], object.wavelike["wavelength"]
        ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
            result.fluxlike["flux"] *= object.fluxlike["flux"]
        else:
            print("Objects do not share wavelength/time axes")
            return

    except AttributeError:

        if type(object) == float or type(object) == int:  # Float or integer.
            result.fluxlike["flux"] *= object

        elif self.shape == object.shape:  # fluxlike
            result.fluxlike["flux"] *= object

        elif len(object) == len(
            self.wavelike["wavelength"]
        ):  # Add along wavelength axis
            result.fluxlike["flux"] *= np.transpose([object] * self.shape[1])
            if self.nwave == self.ntime:
                raise RuntimeError(
                    f"{self} has same number of wavelengths and times; we can't tell which is which."
                )

        elif len(object) == len(self.timelike["time"]):  # Add along time axis.
            result.fluxlike["flux"] *= np.tile(object, (self.shape[0], 1))

        else:
            print("Invalid shape in multiplication: " + str(object.shape))
            return
    return result


def __truediv__(self, object):
    """
    Divide the flux of a rainbow and an input array (or another rainbow)
    and output in a new rainbow object.
    Currently only supports flux division.

    Parameters
    ----------
    object : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow object with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow() with the same parameters as self but with divided
    flux.
    """

    # Create new Rainbow() to store results in.
    result = self._create_copy()

    try:  # Object is another rainbow.
        if np.array_equal(
            self.wavelike["wavelength"], object.wavelike["wavelength"]
        ) and np.array_equal(self.timelike["time"], object.timelike["time"]):
            result.fluxlike["flux"] /= object.fluxlike["flux"]
        else:
            print("Objects do not share wavelength/time axes")
            return

    except AttributeError:

        if type(object) == float or type(object) == int:  # float
            result.fluxlike["flux"] /= object

        elif self.shape == object.shape:  # fluxlike array
            result.fluxlike["flux"] /= object

        elif len(object) == len(
            self.wavelike["wavelength"]
        ):  # Add along wavelength axis
            result.fluxlike["flux"] /= np.transpose([object] * self.shape[1])
            if self.nwave == self.ntime:
                raise RuntimeError(
                    f"{self} has same number of wavelengths and times; we can't tell which is which."
                )
        elif len(object) == len(self.timelike["time"]):  # Add along time axis.
            result.fluxlike["flux"] /= np.tile(object, (self.shape[0], 1))

        else:
            print("Invalid shape in division: " + str(object.shape))
            return
    return result


def __eq__(self, other):
    """
    Test whether (self) == (other).
    """
    # start by assuming the Rainbows are identical
    same = True

    # loop through the core dictionaries
    for d in self._core_dictionaries:

        # pull out each core dictionary from both
        d1, d2 = vars(self)[d], vars(other)[d]

        # loop through elements of each dictionary
        for k in d1:

            # test that all elements match for both
            same *= np.all(d1[k] == d2.get(k, None))

    return same
