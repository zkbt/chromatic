# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_eureka_LCData", "from_eureka_S4"]


def from_eureka_LCData(rainbow, filepath):
    """
    Populate a Rainbow from Eureka! S4 LC outputs,
    which are wavelength-binned light curves.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.

    filepath : str
        The path to the file to load.

    """

    try:
        import astraeus.xarrayIO as xrio
    except ModuleNotFoundError:
        cheerfully_suggest(
            f"""
        You are trying import `astraeus`, which is needed for
        reading and writing `Eureka!` pipeline files. We're
        having trouble importing it, so you could please
        confirm that it is installed into you current
        environment by running the following command?

        pip install git+https://github.com/kevin218/Astraeus.git

        Once you've installed it, please restart whatever
        code you were trying to run, and hopefully it
        will work! Thanks for your patience!
        """
        )
        return None

    # use Astraeus to read in the file
    dataset = xrio.readXR(filepath, verbose=False)

    keys_used = []

    # populate a 1D array of wavelengths (with astropy units of length)
    k = "wavelength"
    w_without_unit = dataset[k].data
    w_unit = u.Unit(dataset[k].attrs["wave_units"])
    rainbow.wavelike["wavelength"] = w_without_unit * w_unit
    keys_used.append(k)

    for k_theirs, k_ours in zip(
        ["wave_hi", "wave_low"], ["wavelength_lower", "wavelength_lower"]
    ):
        w_without_unit = dataset[k_theirs].data
        rainbow.wavelike[k_ours] = w_without_unit * w_unit
        keys_used.append(k_theirs)

    # populate a 1D array of times (with astropy units of time)
    k = "time"
    t_without_unit = dataset[k].data
    if dataset[k].attrs["time_units"] == "BJD_TDB":
        t_unit = u.day
    else:
        t_unit = u.Unit(dataset[t_key].attrs["time_units"])
    rainbow.timelike["time"] = t_without_unit * t_unit
    keys_used.append(k)

    # populate a 2D (row = wavelength, col = time) array of fluxes
    for k in dataset.data_vars.keys():
        shape = dataset[k].shape
        if shape == rainbow.shape:
            q_without_unit = dataset[k].data
            for k_unit in ["flux_units"]:
                try:
                    q_unit = u.Unit(dataset[k_unit])
                    continue
                except (KeyError, ValueError):
                    q_unit = u.Unit()
            rainbow.fluxlike[k] = q_without_unit * q_unit
            keys_used.append(k)
        test_shapes = (rainbow.shape, rainbow.wavelength.shape, rainbow.time.shape)

    rainbow.fluxlike["flux"] = rainbow.data * 1
    rainbow.fluxlike["uncertainty"] = rainbow.err * 1

    keys_used.extend(["wave_mid", "wave_err"])
    for k in dataset.keys():
        if (shape in test_shapes) or (shape[::-1] in test_shapes):
            if k not in keys_used:
                cheerfully_suggest(
                    f"""
                While reading a Eureka S4 LCData dataset, the key '{k}'
                has not been stored as `fluxlike`, `wavelike`, or `timelike`,
                yet its shape {shape} feels like it might match one
                of those. If you'd like to read in that quantity,
                please submit an Issue to the `chromatic` github
                to make sure it gets added!
                """
                )


from_eureka_S4 = from_eureka_LCData
