from ...imports import *


def _is_probably_rainbow(x):
    """
    A quick check if an object is probably a Rainbow
    (without having to do a recursive import to compare
    directly to the class definition).
    """
    return "Rainbow" in x.__class__.__name__


def _do_rainbows_match(a, b):
    """
    A check that two Rainbows have matching wavelengths and times.
    If either of the others aren't Rainbows, an AssertionError
    will be raised.

    Parameters
    ----------
    a : Rainbow
        One Rainbow other.
    b : Rainbow
        Another Rainbow other.

    Returns
    -------
    match : bool
        Do their wavelengths and times match?
    """

    return np.array_equal(a.wavelength, b.wavelength) and np.array_equal(a.time, b.time)


def _raise_ambiguous_shape_error(self, x):
    """
    Raise an error if the shape of a Rainbow is ambiguous.
    """
    if x.nwave == x.ntime:
        raise RuntimeError(
            f"{self} has same number of wavelengths and times; we can't tell which is which."
        )


def _broadcast_to_fluxlike(self, x):
    """
    Convert a scalar, wavelike array, timelike array, or
    fluxlike array into something that can be broadcast into
    a fluxlike shape for math.

    Parameters
    ----------
    x : int, float, bool, array, u.Quantity, ...
        The input.

    Returns
    -------
    x[:,:] : array
        The input in a format that can happily do
        math with a 2D fluxlike array.
    """
    if np.shape(x) == (self.nwave,):
        self._raise_ambiguous_shape_error(self)
        return x[:, np.newaxis]
    elif np.shape(x) == (self.ntime,):
        self._raise_ambiguous_shape_error(self)
        return x[np.newaxis, :]
    elif np.shape(x) in [(), (1,), self.shape]:
        return x
    else:
        raise ValueError(
            f"The shapes {np.shape(x)} and {self.shape} cannot be cast together."
        )


def _apply_operation(self, other, operation, dzdx="1", dzdy="1"):
    """
    Apply a mathematical operation between `self` (x) and `other` (y)
    that returns a new Rainbow where the `flux` (and maybe `model`)
    are the result of the operation (+, -, *, /). Uncertainties
    will be propagated using the string expressions for the required
    derivatives.

    Parameters
    ----------
    self : Rainbow
        The current Rainbow object.
    other : Rainbow, array, float, int
        There are multiple options for what this could be:
            1) scalar float or int
            2) 1D array with same length as wavelength axis
            3) 1D array with same length as time axis
            4) 2D array with same shape as rainbow flux
            5) Rainbow other with same dimensions as self.
    operation : function
        The function that will be applied to combine `self` and `other`.
        Likely options are `np.add`, `np.subtract`, `np.multiply`, and
        `np.true_divide`.
    dzdx : string
        For propagating uncertainties, if z(x,y) is the result,
        what is `dz/dx`? This must be supplied as a string that
        uses the variables `x` (= `self`) and `y` (= `other`).
    dzdy : string
        For propagating uncertainties, if z(x,y) is the result,
        what is `dz/dy`? This must be supplied as a string that
        uses the variables `x` (= `self`) and `y` (= `other`).
    """
    # create new Rainbow() to store results in.
    result = self._create_copy()

    # other object is not Rainbow
    if _is_probably_rainbow(other):
        if _do_rainbows_match(self, other):
            for k in self._keys_that_respond_to_math:
                result.fluxlike[k] = operation(self.fluxlike[k], other.fluxlike[k])
        else:
            raise ValueError(
                f"The two Rainbow objects {self} and {other} don't share wavelength/time axes."
            )
        sigma_y = other.uncertainty
        y = other.get("model")
        if y is None:
            y = other.flux

    # other object not Rainbow
    else:
        y = self._broadcast_to_fluxlike(other)
        sigma_y = 0
        for k in self._keys_that_respond_to_math:
            result.fluxlike[k] = operation(self.fluxlike[k], y)

    sigma_x = self.uncertainty
    x = self.get("model")
    if x is None:
        x = self.flux

    # If z = operation(x,y), then to propagate errors we need
    # to use the derivatives of z with respect to x (dz/dx) and y (dz/dy):
    # print(f"mean(x) = {np.mean(x)}")
    # print(f"mean(sigma_x) = {np.mean(sigma_x)}")
    # print(f"dzdx = {dzdx}")
    # print(f"mean(y) = {np.mean(y)}")
    # print(f"mean(sigma_y) = {np.mean(sigma_y)}")
    # print(f"dzdy = {dzdy}")

    variance = sigma_x**2 * eval(dzdx) ** 2 + sigma_y**2 * eval(dzdy) ** 2
    result.fluxlike["uncertainty"] = np.sqrt(variance)

    return result


def __add__(self, other):
    """
    Add the flux of a rainbow and an input array (or another rainbow)
    and output in a new rainbow other.

    Parameters
    ----------
    other : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow other with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow
        A new `Rainbow` with the mathematical operation applied.
    """

    # create the history entry
    h = self._create_history_entry("+", locals())

    # calculate a new Rainbow using the operation and error propagation
    result = self._apply_operation(other, operation=np.add, dzdx="1", dzdy="1")

    # append the history entry to the new Rainbow
    result._record_history_entry(h)

    return result


def __sub__(self, other):
    """
    Subtract the flux of a rainbow from an input array (or another rainbow)
    and output in a new rainbow other.

    Parameters
    ----------
    other : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow other with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow
        A new `Rainbow` with the mathematical operation applied.
    """
    # create the history entry
    h = self._create_history_entry("-", locals())

    # calculate a new Rainbow using the operation and error propagation
    result = self._apply_operation(other, operation=np.subtract, dzdx="1", dzdy="1")

    # append the history entry to the new Rainbow
    result._record_history_entry(h)

    return result


def __mul__(self, other):
    """
    Multiply the flux of a rainbow and an input array (or another rainbow)
    and output in a new rainbow other.

    Parameters
    ----------
    other : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow other with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow
        A new `Rainbow` with the mathematical operation applied.
    """

    # create the history entry
    h = self._create_history_entry("*", locals())

    # calculate a new Rainbow using the operation and error propagation
    result = self._apply_operation(other, operation=np.multiply, dzdx="y", dzdy="x")

    # append the history entry to the new Rainbow
    result._record_history_entry(h)

    return result


def __truediv__(self, other):
    """
    Divide the flux of a rainbow and an input array (or another rainbow)
    and output in a new rainbow other.

    Parameters
    ----------
    other : Array or float.
        Multiple options:
        1) float
        2) 1D array with same length as wavelength axis
        3) 1D array with same length as time axis
        4) 2D array with same shape as rainbow flux
        5) Rainbow other with same dimensions as self.

    Returns
    ----------
    rainbow : Rainbow
        A new `Rainbow` with the mathematical operation applied.
    """
    # create the history entry
    h = self._create_history_entry("/", locals())

    # calculate a new Rainbow using the operation and error propagation
    result = self._apply_operation(
        other, operation=np.true_divide, dzdx="1/y", dzdy="-x/y**2"
    )

    # append the history entry to the new Rainbow
    result._record_history_entry(h)

    return result


def __eq__(self, other):
    """
    Test whether `self == other`.

    This compares the wavelike, timelike, and fluxlike arrays
    for exact matches. It skips entirely over the metadata.

    Parameters
    ----------
    other : Rainbow
        Another `Rainbow` to compare to.

    Returns
    -------
    equal : bool
        Are they (effectively) equivalent?
    """
    # start by assuming the Rainbows are identical
    same = True

    for a, b in zip([self, other], [other, self]):

        # loop through the core dictionaries
        for d in a._core_dictionaries:
            if d == "metadata":
                continue
            # pull out each core dictionary from both
            d1, d2 = vars(self)[d], vars(b)[d]
            same *= set(d1.keys()) == set(d2.keys())

            # loop through elements of each dictionary
            for k in d1:

                # ignore different histories (e.g. new vs loaded)
                if k != "history":

                    # test that all elements match for both
                    if d == "fluxlike":
                        same *= np.all(
                            np.isclose(
                                a.get(k)[a.ok.astype(bool)], b.get(k)[a.ok.astype(bool)]
                            )
                        )
                    else:
                        same *= np.all(np.isclose(a.get(k), b.get(k)))

    return bool(same)


def diff(self, other):
    """
    Test whether `self == other`, and print the differences.

    This compares the wavelike, timelike, and fluxlike arrays
    for exact matches. It skips entirely over the metadata.
    The `diff` function is the same as `__eq__`, but a little
    more verbose, just to serve as a helpful debugging tool.

    Parameters
    ----------
    other : Rainbow
        Another `Rainbow` to compare to.

    Returns
    -------
    equal : bool
        Are they (effectively) equivalent?
    """
    # start by assuming the Rainbows are identical
    same = True

    for a, b in zip([self, other], [other, self]):

        # loop through the core dictionaries
        for d in a._core_dictionaries:
            if d == "metadata":
                continue
            # pull out each core dictionary from both
            d1, d2 = vars(self)[d], vars(b)[d]
            if set(d1.keys()) != set(d2.keys()):
                differences = list(set(d1.keys()) - set(d2.keys()))
                print(f"{a}.{d} has {differences} and {b} does not")

            # loop through elements of each dictionary
            for k in d1:
                # ignore different histories (e.g. new vs loaded)
                if k != "history":
                    continue
                # test that all elements match for both
                if np.all(a.get(k) == b.get(k)) == False:
                    print(f"{a}.{d}[{k}] != {b}.{d}[{k}]")
