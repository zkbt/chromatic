# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_hst_stis"]


def get_wavelengths_from_x1dints_files(filenames, order=1, end_wavelengths_to_remove=0):
    """
    Wrapper to extract wavelengths from a list of file segments.

    Parameters
    ----------
    filenames : list
        The files to open, hopefully sorted.
    order : int or list or None
        The orders to extract, default=1
    end_wavelengths_to_remove : int
        The number of 'edge' wavelengths to remove from each order, default=0

    """
    with fits.open(filenames[0]) as hdu:
        unit_string = hdu["SCI"].columns["wavelength"].unit
        if unit_string is None:
            unit_string = "micron"
            cheerfully_suggest("No wavelength unit was found; assuming 'micron'.")
        try:
            wavelength_unit = u.Unit(unit_string)
        except:
            # remove plural from unit
            if unit_string[-1] == "s":
                wavelength_unit = u.Unit(unit_string[:-1])

        wave_shape = np.shape(hdu["SCI"].data["wavelength"])

        if order is None:
            n_orders = wave_shape[0]
            wavelengths_per_order = wave_shape[1]
            wave, not_ok = [], []
            for n in range(n_orders):
                w = hdu["SCI"].data["wavelength"][n]
                wave.extend(w)
                not_ok.extend(
                    1
                    * (
                        (np.arange(wave_shape[1]) >= end_wavelengths_to_remove)
                        & (
                            np.arange(wave_shape[1])
                            < wave_shape[1] - end_wavelengths_to_remove
                        )
                    )
                )
        elif type(order) == int:
            wave = hdu["SCI"].data["wavelength"][order]
            not_ok = 1 * (
                (np.arange(wave_shape[1]) >= end_wavelengths_to_remove)
                & (np.arange(wave_shape[1]) < wave_shape[1] - end_wavelengths_to_remove)
            )
        else:
            wave, not_ok = [], []
            for n, order_n in enumerate(order):
                w = hdu["SCI"].data["wavelength"][order_n]
                wave.extend(w)
                not_ok.extend(
                    1
                    * (
                        (np.arange(wave_shape[1]) >= end_wavelengths_to_remove)
                        & (
                            np.arange(wave_shape[1])
                            < wave_shape[1] - end_wavelengths_to_remove
                        )
                    )
                )

    return {"wavelength": wave * wavelength_unit * 1}, not_ok  # , "order":orders}


def get_times_from_x1dints_files(filenames):
    """
    Wrapper to extract or estimate times from a list of file segments.

    Parameters
    ----------
    filenames : list
        The files to open, hopefully sorted.

    """
    # first, try to get times from the `TEXPSTRT` and `TEXPEND` headers
    try:
        timelike = {}
        for f in filenames:
            with fits.open(f) as hdu:
                t_start = hdu["PRIMARY"].header["TEXPSTRT"]
                t_end = hdu["PRIMARY"].header["TEXPEND"]
                if "mjd" in timelike:
                    timelike["mjd"].extend(list([0.5 * (t_start + t_end)]))
                else:
                    timelike["mjd"] = list([0.5 * (t_start + t_end)])

        for c in timelike:
            if "JD" in c:
                timelike[c] = timelike[c] * u.day
            else:
                timelike[c] = np.asarray(timelike[c])

        #         timelike['time'] = timelike['mjd']

        return timelike

    except KeyError:
        cheerfully_suggest(
            f"""
        No `EXPSTART/EXPEND` extension was found in the first file
        {os.path.basename(filenames[0])}
        """
        )


def from_hst_stis(
    rainbow, filepath, order=None, different_visits=True, end_wavelengths_to_remove=0
):
    """
    Populate a Rainbow from a file in the HST format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated. This function is meant
        to be called during the initialization of a Rainbow
        object. You can/should assume that the `rainbow` object
        being passed will already contain the four core
        dictionaries `timelike`, `wavelike`, `fluxlike`, and
        `metadata`. This function should, at a minimum, add
        the following items
            + `self.timelike['time']`
            + `self.wavelike['wavelength']`
            + `self.fluxlike['flux']`
        and optionally, additional entries like
            + `self.metadata['some-useful-parameter']`
            + `self.fluxlike['uncertainty']`
            + `self.fluxlike['ok']`

    filepath : str
        The path to the file to load.

    order : None or int or list
        The order(s) of wavelengths to load. If None assume
        we want to import all orders into one Rainbow

    different_visits : Boolean
        Whether we want to store the headers/metadata from each
        time rather than just the first.
    """

    # create a list of filenames (which might be just 1)
    filenames = expand_filenames(filepath)

    # loop over file (each one a segment)
    for i_file, f in enumerate(filenames):
        # open this fits file
        with fits.open(f) as hdu:

            if different_visits:
                # if we have different visits assume they will have different headers
                if "header" in rainbow.metadata.keys():
                    rainbow.metadata["header"][i_file] = hdu["PRIMARY"].header
                else:
                    rainbow.metadata["header"] = {i_file: hdu["PRIMARY"].header}

            if i_file == 0:

                if not different_visits:
                    # just use one header for all data
                    rainbow.metadata["header"] = hdu["PRIMARY"].header

                wave_shape = np.shape(hdu["SCI"].data["wavelength"])

                # find the number of orders
                if len(wave_shape) > 1:
                    n_orders = wave_shape[0]
                    wavelengths_per_order = wave_shape[1]
                else:
                    n_orders = 1
                    wavelengths_per_order = wave_shape[0]

                    if order is None:
                        if n_orders > 1:
                            cheerfully_suggest(
                                f"""
                                This file contains data for {n_orders} spectrosopic orders. Because no
                                `order=` keyword was supplied, we're defaulting to importing all orders. You can
                                hide this warning by expliciting stating which order you want to load.
                                For this file, the options include {np.arange(n_orders) + 1}.
                                """
                            )
                #                 assert (order >= 1) and (order <= n_orders)

                wavelike, not_ok = get_wavelengths_from_x1dints_files(
                    filenames,
                    order=order,
                    end_wavelengths_to_remove=end_wavelengths_to_remove,
                )
                rainbow.wavelike.update(**wavelike)
                rainbow.wavelike["ok"] = np.array(not_ok)

                timelike = get_times_from_x1dints_files(filenames)
                rainbow.timelike.update(**timelike)
                astropy_times = Time(timelike["mjd"], format="mjd", scale="tdb")
                rainbow.set_times_from_astropy(astropy_times, is_barycentric=True)

                column_units = {}
                for column in hdu["SCI"].columns:

                    # get a lower case name for the unit
                    c = column.name.lower()

                    # set up an empty Quantity array, with units if known
                    this_quantity = u.Quantity(
                        np.nan * np.ones((rainbow.nwave, rainbow.ntime))
                    )
                    if column.unit is not None:
                        try:
                            column_units[c] = u.Unit(column.unit)
                        except:
                            try:
                                # remove plural from unit
                                if column.unit[-1] == "s":
                                    column_units[c] = u.Unit(column.unit[:-1])
                            except:
                                column_units[c] = 1

                        this_quantity *= column_units[c]
                    else:
                        column_units[c] = 1

                    # populate the fluxlike dictionary with the empty array
                    rainbow.fluxlike[c] = this_quantity * 1

            # loop through all the columns in the data extension
            for column in hdu["SCI"].columns:

                # get a lower case name for the unit
                c = column.name.lower()
                current_time_index = i_file

                if order is None:
                    for n in range(n_orders):
                        # store in the appropriate column of fluxlike array
                        sci_data = hdu["SCI"].data[c][n]
                        rainbow.fluxlike[c][
                            n
                            * int(rainbow.nwave / n_orders) : (n + 1)
                            * int(rainbow.nwave / n_orders),
                            current_time_index,
                        ] = (
                            sci_data * column_units[c] * 1
                        )
                elif type(order) == int:
                    # store in the appropriate column of fluxlike array
                    sci_data = hdu["SCI"].data[c][order - 1]
                    rainbow.fluxlike[c][:, current_time_index] = (
                        sci_data * column_units[c] * 1
                    )
                else:
                    n_orders = len(order)
                    for n, order_n in enumerate(order):
                        # store in the appropriate column of fluxlike array
                        sci_data = hdu["SCI"].data[c][order_n - 1]
                        rainbow.fluxlike[c][
                            n
                            * int(rainbow.nwave / n_orders) : (n + 1)
                            * int(rainbow.nwave / n_orders),
                            current_time_index,
                        ] = (
                            sci_data * column_units[c] * 1
                        )

    if "uncertainty" not in rainbow.fluxlike.keys():
        if "error" in rainbow.fluxlike.keys():
            rainbow.fluxlike["uncertainty"] = rainbow.fluxlike["error"]
            rainbow.fluxlike.pop("error")

    rainbow.fluxlike["ok"] = rainbow.fluxlike["uncertainty"] > 0
