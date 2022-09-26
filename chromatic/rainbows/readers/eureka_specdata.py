# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_eureka_SpecData", "from_eureka_S3"]


def from_eureka_SpecData(rainbow, filepath, optimal=True):
    """
    Populate a Rainbow from Eureka! S3/S4 SpecData outputs,
    which are time-series extracted spectra.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.

    filepath : str
        The path to the file to load.

    optimal : bool
        Should the optimal extraction be used as the flux?
        If True, (flux, uncertainty) = (optspec, opterr)
        If False, (flux, uncertainty) = (stdspec, sqrt(stdvar))
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
    k = "wave_1d"
    w_without_unit = dataset[k].data
    w_unit = u.Unit(dataset[k].attrs["wave_units"])
    rainbow.wavelike["wavelength"] = w_without_unit * w_unit
    keys_used.append(k)

    # populate a 1D array of times (with astropy units of time)
    k = "time"
    if dataset[k].attrs["time_units"] == "BJD_TDB":
        astropy_times = Time(dataset[k].data, format="mjd", scale="tdb")
        rainbow.set_times_from_astropy(astropy_times, is_barycentric=True)
    else:
        t_unit = u.Unit(dataset[t_key].attrs["time_units"])
        t_without_unit = dataset[k].data
        rainbow.timelike["time"] = t_without_unit * t_unit
    keys_used.append(k)

    # populate a 2D (row = wavelength, col = time) array of fluxes
    for k in dataset.data_vars.keys():
        shape = dataset[k].shape
        if shape[::-1] == rainbow.shape:
            q_without_unit = dataset[k].data.T
            for k_unit in ["flux_units"]:
                try:
                    q_unit = u.Unit(dataset[k_unit])
                    continue
                except (KeyError, ValueError):
                    q_unit = u.Unit()
            rainbow.fluxlike[k] = q_without_unit * q_unit
            keys_used.append(k)
        test_shapes = (rainbow.shape, rainbow.wavelength.shape, rainbow.time.shape)

    if optimal:
        rainbow.fluxlike["flux"] = rainbow.optspec * 1
        rainbow.fluxlike["uncertainty"] = rainbow.opterr * 1
    else:
        rainbow.fluxlike["flux"] = rainbow.stdspec * 1
        rainbow.fluxlike["uncertainty"] = np.sqrt(rainbow.stdvar) * 1

    for k in dataset.keys():
        if k not in keys_used:
            shape = dataset[k].shape
            values = dataset[k].data
            try:
                if shape in test_shapes:
                    rainbow._put_array_in_right_dictionary(k, values)
                elif shape[::-1] in test_shapes:
                    rainbow._put_array_in_right_dictionary(k, values.T)
            except ValueError:
                cheerfully_suggest(
                    f"""
                While reading a Eureka S3 SpecData dataset, the key '{k}'
                has not been stored as `fluxlike`, `wavelike`, or `timelike`,
                but its shape of {shape} feels like maybe it should fit.
                If you'd like to read in that quantity, please submit an
                Issue to the `chromatic` github to make sure it gets added!
                """
                )


from_eureka_S3 = from_eureka_SpecData
