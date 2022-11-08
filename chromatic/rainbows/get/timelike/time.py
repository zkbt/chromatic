from ....imports import *


def get_times_as_astropy(self, time=None, format=None, scale=None, is_barycentric=None):
    """
    Convert times from a `Rainbow` into an astropy `Time` object.

    Parameters
    ----------
    time : Quantity, optional
        The time-like Quantity to be converted.
        If None (default), convert the time values in `self.time`
        If another time-like Quantity, convert those values.
    format : str, optional
        The time format to supply to astropy.time.Time.
        If None (default), format will be pulled from
        `self.metadata['time_details']['format']`
    scale : str, optional
        The time scale to supply to astropy.time.Time.
        If None (default), scale will be pulled from
        `self.metadata['time_details']['scale']`
    is_barycentric : bool, optional
        Are the times already measured relative to the
        Solar System barycenter? This is mostly for warning
        the user that it's not.
        If `None` (default), `is_barycentric` will be pulled from
        `self.metadata['time_details']['is_barycentric']`

    Returns
    -------
    astropy_time : Time
        The times as an astropy `Time` object.
    """

    # take times from self or from the keyword
    if time is None:
        time = self.time

    # give a format warning
    format = format or self.get("time_format")
    if format is None:
        cheerfully_suggest(
            f"""
        `.metadata['time_details']['format']` is not set,
        nor was a `format=` keyword argument provided.

        Since `.time` is already an astropy Quantity,
        this is likely a question of whether the format is
        'jd' or 'mjd' (= 'jd' - 2400000.5) or something else.
        If you get this wrong, you might be lost in time!

        For more about astropy.Time formats, please see:
        https://docs.astropy.org/en/stable/time/index.html#time-format
        """
        )

    # give a scale warning
    scale = scale or self.get("time_scale")
    if scale is None:
        now = Time.now()
        differences_string = ""
        for s in now.SCALES:
            dt = ((getattr(now, s).jd - now.tdb.jd) * u.day).to(u.second)
            differences_string += f"{s:>15} - tdb = {dt:10.6f}\n"
        cheerfully_suggest(
            f"""
        .metadata['time_details']['scale'] is not set,
        nor was a `scale=` keyword argument provided.

        The main question is whether the time scale is 'tdb'
        (Barycentric Dynamical Time) or something close to it,
        or 'utc' or something close to it. The differences
        between these options, at {now.utc.iso} (UTC), are:
        \n{differences_string}
        If you get this wrong, you might be lost in time!

        For more about astropy.Time scales, please see:
        https://docs.astropy.org/en/stable/time/index.html#time-scale
        """
        )

    # give some barycenter warnings
    is_barycentric = is_barycentric or self.get("time_is_barycentric")
    if is_barycentric == True and "ut" in scale.lower():
        cheerfully_suggest(
            f"""
        barycentic={is_barycentric} and scale={scale}
        It's a deeply weird combination to have a barycentric
        time measured at the Solar System barycentric but in
        Earth's leap-second-based UTC system. Please consider
        checking your time details.
        """
        )
    if is_barycentric != True:
        cheerfully_suggest(
            f"""
        The returned time is not known to be measured relative
        to the Solar System barycenter. It's probably therefore
        measured from Earth or the position of your telescope,
        but please be warned that the timing of very distant event
        (like exoplanet transits) might be off by up to about
        8 minutes (= the light travel time between Earth + Sun).
        """
        )

    # generate astropy Time array
    astropy_time = Time(self.time, format=format, scale=scale)

    # do a check that the values aren't really weird
    if (astropy_time.min().decimalyear < 1000) or (
        astropy_time.max().decimalyear > 3000
    ):
        cheerfully_suggest(
            f"""
        The times, which span
        jd={astropy_time.min().jd} to jd={astropy_time.max().jd}
        don't seem likely to be within the range of modern astronomical
        observations. Please consider double checking your time values
        and/or (format='{format}', scale='{scale}').
        """
        )

    return astropy_time


def set_times_from_astropy(self, astropy_time, is_barycentric=None):
    """
    Set the times for this `Rainbow` from an astropy `Time` object.

    Parameters
    ----------
    astropy_time : Time
        The times as an astropy `Time` object.
    is_barycentric : bool, optional
        Are the times already measured relative to the
        Solar System barycenter? This is mostly for warning
        the user that it's not. Options are True, False,
        None (= don't know).

    Returns
    -------
    time : Quantity
        An astropy Quantity with units of time,
        expressing the Time as julian day.
        In addition to this returned variable,
        the function sets the following internal
        variables:
        ```
        self.time # (= the astropy Quantity of times)
        self.metadata['time_format'] # (= the format to convert back to Time)
        self.metadata['time_scale'] # (= the scale to convert back to Time)
        self.metadata['time_is_barycentric'] # (= is it barycentric?)
        ```
    """

    # set the formats
    format = "jd"
    unit = u.day
    scale = "tdb"

    # store the necessary values
    self.timelike["time"] = getattr(getattr(astropy_time, scale), format) * unit
    self.metadata["time_format"] = format
    self.metadata["time_scale"] = scale
    self.metadata["time_is_barycentric"] = is_barycentric

    # do some accounting to sync everything together
    self._guess_tscale()
    self._make_sure_time_edges_are_defined()
    return self.time
