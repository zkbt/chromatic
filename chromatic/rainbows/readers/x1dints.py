"""
Define a reader for STScI pipeline x1dints.fits files.
"""
from ...imports import *

__all__ = ["from_x1dints"]


def get_times_from_x1dints_files(filenames):
    """
    Wrapper to extract or estimate times from a list of file segments.

    Parameters
    ----------
    filenames : list
        The files to open, hopefully sorted.

    """
    try:
        timelike = {}
        for f in filenames:
            with fits.open(f) as hdu:
                table_of_times = hdu["int_times"].data
                for c in table_of_times.columns.names:
                    if c in timelike:
                        timelike[c].extend(list(table_of_times[c]))
                    else:
                        timelike[c] = list(table_of_times[c])

        for c in timelike:
            if "JD" in c:
                timelike[c] = timelike[c] * u.day
            else:
                timelike[c] = np.asarray(timelike[c])

        return timelike

    except KeyError:
        cheerfully_suggest(
            f"""
        No `int_times` extension was found in the first file
        {filenames[0]}
        We're estimating the times from the `sci` extension
        using the TDB-BEG, TDB-END, and EFFINTTM keywords.
        """
        )

    try:
        with fits.open(filenames[0]) as hdu:

            # make up a time grid (this is still a little kludgy; it'd be better to have int_times)
            first_mjd_barycentric_integration_midpoint = (
                hdu["sci"].header["TDB-BEG"] * u.day
                + hdu["primary"].header["EFFINTTM"] / 2 * u.second
            )
            last_mjd_barycentric_integration_midpoint = (
                hdu["sci"].header["TDB-END"] * u.day
                - hdu["primary"].header["EFFINTTM"] / 2 * u.second
            )

        n_integrations_total = hdu["PRIMARY"].header["NINTS"]

        mjd_barycentric_integration_midpoints = np.linspace(
            first_mjd_barycentric_integration_midpoint,
            last_mjd_barycentric_integration_midpoint,
            n_integrations_total,
        )
        cheerfully_suggest(
            f"""
        Times were set by linearly interpolating between the exposure
        start and end points. It's very possible these times are off
        by at least a few seconds and possibly up to the duration
        of one integration (= {hdu["primary"].header["EFFINTTM"]}s)."""
        )
        return {"int_mid_BJD_TDB": mjd_barycentric_integration_midpoints}
    except KeyError:
        raise ValueError(
            f"""
        Blerg! We can't seem to determine the integration times
        from the `INT_TIMES` extension or even estimate them from
        the `SCI` extension.
        """
        )


def get_wavelengths_from_x1dints_files(
    filenames, order=1, non_spectrum_extensions=["PRIMARY", "INT_TIMES"]
):

    with fits.open(filenames[0]) as hdu:

        #
        e = (len(non_spectrum_extensions) - 1) + order - 1
        unit_string = hdu[e].columns["wavelength"].unit
        if unit_string is None:
            unit_string = "micron"
            cheerfully_suggest("No wavelength unit was found; assuming 'micron'.")
        wavelength_unit = u.Unit(unit_string)
        return {"wavelength": hdu[e].data["wavelength"] * wavelength_unit * 1}


def from_x1dints(rainbow, filepath, order=1, **kw):
    """
    Populate a Rainbow from an STScI pipeline x1dints file,
    or a group of x1dints files for multiple segments.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated. (This is intended
        to be used as a class method of Rainbow or
        of a class derived from Rainbow, as a way of
        initializing an object from files.)
    filepath : str
        The path to the file, or a glob-friendly string
        pointing to a group of files that should all be
        loaded together (for example, if an exposure was
        split into multiple segments).
    """

    # create a list of filenames (which might be just 1)
    filenames = expand_filenames(filepath)

    # loop over file (each one a segment)
    for i_file, f in enumerate(filenames):

        # open this fits file
        with fits.open(f) as hdu:

            # if this is the first one, populate the shared stuff
            if i_file == 0:

                # grab the entire primary header
                rainbow.metadata["header"] = hdu["PRIMARY"].header

                if hdu["PRIMARY"].header["INSTRUME"] == "NIRISS":
                    n_orders = 3
                    cheerfully_suggest(
                        f"""
                    Loading NIRISS spectroscopic `order={order}``. Three orders are available,
                    and you can set which (1,2,3) you want to read with the `order=` option.
                    """
                    )
                else:
                    n_orders = 1
                assert (order >= 1) and (order <= n_orders)

                non_spectrum_extensions = [
                    x for x in ["PRIMARY", "SCI", "INT_TIMES", "ASDF"] if x in hdu
                ]

            # set the index to the start of this segment
            integration_counter = hdu["PRIMARY"].header["INTSTART"]
            n_integrations_in_this_segment = (
                hdu["PRIMARY"].header["INTEND"] - hdu["PRIMARY"].header["INTSTART"] + 1
            )

            # make sure sizes match, ignoring PRIMARY, SCI, ASDF
            assert n_integrations_in_this_segment * n_orders == (
                len(hdu) - len(non_spectrum_extensions)
            )

            if i_file == 0:

                wavelike = get_wavelengths_from_x1dints_files(
                    filenames,
                    order=order,
                    non_spectrum_extensions=non_spectrum_extensions,
                )
                rainbow.wavelike.update(**wavelike)

                timelike = get_times_from_x1dints_files(filenames)
                rainbow.timelike.update(**timelike)
                astropy_times = Time(
                    timelike["int_mid_BJD_TDB"], format="mjd", scale="tdb"
                )
                rainbow.set_times_from_astropy(astropy_times, is_barycentric=True)

            # loop through the integrations in this segment
            for i in tqdm(np.arange(n_integrations_in_this_segment), leave=False):
                e = 2 + i * n_orders + (order - 1)

                # do this stuff only on the very first time through
                if i_file == 0 and i == 0:

                    # set up the fluxlike quantities
                    column_units = {}
                    for column in hdu[e].columns:

                        # get a lower case name for the unit
                        c = column.name.lower()

                        # set up an empty Quantity array, with units if known
                        this_quantity = u.Quantity(
                            np.nan * np.ones((rainbow.nwave, rainbow.ntime))
                        )
                        if column.unit is not None:
                            column_units[c] = u.Unit(column.unit)
                            this_quantity *= column_units[c]
                        else:
                            column_units[c] = 1

                        # populate the fluxlike dictionary with the empty array
                        rainbow.fluxlike[c] = this_quantity

                for column in hdu[e].columns:

                    # get a lower case name for the unit
                    c = column.name.lower()
                    rainbow.fluxlike[c][:, integration_counter - 1] = (
                        hdu[e].data[c] * column_units[c]
                    )

                # increment the running integration total
                integration_counter += 1

    n_filled_times = np.sum(np.any(np.isfinite(rainbow.flux), rainbow.waveaxis))
    if n_filled_times != rainbow.ntime:
        cheerfully_suggest(
            f"""
        The x1dints header(s) indicate there should be {rainbow.ntime} integrations,
        but only {n_filled_times} columns of the flux array were populated. Are you
        perhaps missing some segment files?
        """
        )

    # try to pull in the errors
    try:
        rainbow.fluxlike["uncertainty"] = rainbow.fluxlike["flux_error"] * 1
    except KeyError:
        rainbow.fluxlike["uncertainty"] = rainbow.fluxlike["error"] * 1

    if rainbow.uncertainty is None:
        message = f"""
        Hmmm...it's not clear which column corresponds to the
        flux uncertainties for this Rainbow object. The
        available `fluxlike` columns are:
            {rainbow.fluxlike.keys()}
        A long-term solution might be to fix the `from_x1dints`
        reader, but a short-term solution would be to pick one
        of the columns listed above and say something like

        x.fluxlike['uncertainty'] = x.fluxlike['some-other-relevant-error-column']

        where `x` is the Rainbow you just created.
        """
        cheerfully_suggest(message)
