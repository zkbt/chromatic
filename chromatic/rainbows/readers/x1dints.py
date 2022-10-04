"""
Define a reader for STScI pipeline x1dints.fits files.
"""
from ...imports import *

__all__ = ["from_x1dints"]


def guess_jwst_pipeline_stage(hdu):
    """
    Guess whether this file is a Stage 2 or Stage 3 pipeline output.

    Parameters
    ----------
    hdu : FITS Header Data Unit
        The (first) loaded FITS file of a dataset.

    Return
    ------
    stage : int
        2 (= before outlier rejection)
        3 (= after outlier rejection, which still seems a little buggy)
    """
    # according to...
    # https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
    # outputs from stages 0-2 shoud look like
    # jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>(-<”seg”NNN>)_<detector>_<prodType>.fits
    # outputs from stage 3 should look like
    # jw<ppppp>-<AC_ID>_[<”t”TargID | “s”SourceID>](-<”epoch”X>)_<instr>_<optElements>(-<subarray>)_<prodType>(-<ACT_ID>).fits
    filename = hdu["PRIMARY"].header["FILENAME"]
    basename = os.path.basename(filename)
    if "x1dints" not in basename:
        cheerfully_suggest(
            f"""
        The file you're trying to load does not appear to be a 'x1dints' file
        produced by the STScI `jwst` pipeline. Its filename (in the FITS header) is
        '{basename}'
        which does not contain the string 'x1dints'. Please use a different
        chromatic reader format beside 'x1dints' to try to read this file.
        """
        )

    # outlier rejection happens only in stage 3
    if "S_OUTLIR" in hdu["PRIMARY"].header:
        return 3
    else:
        return 2


def get_times_from_x1dints_files(filenames):
    """
    Wrapper to extract or estimate times from a list of file segments.

    Parameters
    ----------
    filenames : list
        The files to open, hopefully sorted.

    """
    # first, try to get times from an `INT_TIMES` extension
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
        {os.path.basename(filenames[0])}
        """
        )

    # second, try to estimate times from a `SCI` or `PRIMARY` extension
    try:
        with fits.open(filenames[0]) as hdu:
            try:
                key = "SCI"
                extension = hdu[key]
            except KeyError:
                key = "PRIMARY"
                extension = hdu[key]

            # figure out the right keyword to get barycentric time
            start_key = [x for x in ["TDB-BEG", "BSTRTIME"] if x in extension.header][0]
            end_key = [x for x in ["TDB-END", "BENDTIME"] if x in extension.header][0]

            # make up a time grid (this is still a little kludgy; it'd be better to have int_times)
            first_mjd_barycentric_integration_midpoint = (
                extension.header[start_key] * u.day
                + hdu["primary"].header["EFFINTTM"] / 2 * u.second
            )
            last_mjd_barycentric_integration_midpoint = (
                extension.header[end_key] * u.day
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
        start and end points, as estimated from the '{key}' extension
        using the '{start_key}', '{end_key}', and 'EFFINTTM' keywords.
        Times may be off by a few seconds and possibly up to the duration
        of one integration (= {hdu["primary"].header["EFFINTTM"]}s)."""
        )
        return {"int_mid_BJD_TDB": mjd_barycentric_integration_midpoints}
    except KeyError:
        raise ValueError(
            f"""
        Blerg! We can't seem to determine the integration times
        from the `INT_TIMES` extension or even estimate them from
        the `{key}` extension. Being unsure how to proceed here,
        it might be a good idea to ask about your dataset at
            https://github.com/zkbt/chromatic/issues
        including the errors and warnings you're finding, along
        with a clear description of what data you're trying on.
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


def from_x1dints(rainbow, filepath, order=None, **kw):
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

                # store the entire primary header in metadata
                rainbow.metadata["header"] = hdu["PRIMARY"].header

                # figure out the number of orders allowed for this dataset
                if hdu["PRIMARY"].header["INSTRUME"] == "NIRISS":
                    n_orders = 3
                else:
                    n_orders = 1

                if order is None:
                    order = 1
                    if n_orders > 1:
                        cheerfully_suggest(
                            f"""
                        This file contains data for {n_orders} spectrosopic orders. Because no
                        `order=` keyword was supplied, we're defaulting to first order. You can
                        hide this warning by expliciting stating which order you want to load.
                        For this file, the options include {np.arange(n_orders) + 1}.
                        """
                        )
                assert (order >= 1) and (order <= n_orders)

                # guess if this file is stage 3 or earlier (= stage 2)
                pipeline_stage = guess_jwst_pipeline_stage(hdu)
                if pipeline_stage == 3:
                    cheerfully_suggest(
                        f"""
                    YIKES! In our testing, we've found that some data products from Stage 3 of the
                    STScI `jwst` pipeline have weird problems, including integration segments that
                    have been stitched together in the wrong temporal order and/or missing time
                    information for the individual integrations. If the `chromatic` reader succeeds
                    in loading your requested data into a `Rainbow` object, you should still be
                    very suspicious of them!

                    A reasonable alternative, if you just want a quick look at the time-series
                    spectra for your dataset, is to try to load the Stage 2 pipeline `x1dints`
                    files. They won't have Stage 3's outlier-rejection applied (of which folks
                    anyway still a little suspicious) but should otherwise be similar. These
                    file(s) may be split into multiple segments, each with a format like
                    `jw02734002001_04101_00001-seg*_nis_x1dints.fits`
                    where the `*` can be used to point to all matching files in a location.
                    """
                    )

                # define a complete list of
                non_spectrum_extensions = [
                    x for x in ["PRIMARY", "SCI", "INT_TIMES", "ASDF"] if x in hdu
                ]
                non_asdf_non_spectrum_extensions = non_spectrum_extensions[:-1]

            # set the index to the start of this segment
            integration_counter = hdu["PRIMARY"].header["INTSTART"]
            n_integrations_in_this_segment = (
                hdu["PRIMARY"].header["INTEND"] - hdu["PRIMARY"].header["INTSTART"] + 1
            )

            # make sure sizes match, ignoring PRIMARY, SCI, ASDF
            N_integration_extensions = len(hdu) - len(non_spectrum_extensions)
            N_expected_integration_extensions = (
                n_integrations_in_this_segment * n_orders
            )

            if pipeline_stage == 2:
                # so far, all real Stage 2 files seem to behave normally
                assert N_expected_integration_extensions == N_integration_extensions
            elif pipeline_stage == 3:
                # some Stage 3 files seem to have errors in the `INTSTART` + `INTEND` keywords
                if N_expected_integration_extensions != N_integration_extensions:
                    cheerfully_suggest(
                        f"""
                    The number of available integration extensions ({N_integration_extensions}) in this file
                    does not match the expected number ({N_expected_integration_extensions}) for this segment.
                    However, this seems to be a Stage 3 `x1dints` file where the
                    `INTSTART` + `INTEND` keywords sometimes behave weirdly, so it's
                    probably OK not to worry about this. We'll proceed by trying
                    to load as many EXTRACT1D extensions as we can."""
                    )
                n_integrations_in_this_segment = int(
                    N_integration_extensions / n_orders
                )

            # get wavelengths and times from the first file
            if i_file == 0:

                # get the wavelengths from the first file (assuming they're constant)
                wavelike = get_wavelengths_from_x1dints_files(
                    filenames,
                    order=order,
                    non_spectrum_extensions=non_spectrum_extensions,
                )
                rainbow.wavelike.update(**wavelike)

                # get the times from the first file
                timelike = get_times_from_x1dints_files(filenames)
                rainbow.timelike.update(**timelike)
                astropy_times = Time(
                    timelike["int_mid_BJD_TDB"], format="mjd", scale="tdb"
                )
                rainbow.set_times_from_astropy(astropy_times, is_barycentric=True)

            # loop through the integrations in this segment
            for i in tqdm(np.arange(n_integrations_in_this_segment), leave=False):
                e = len(non_asdf_non_spectrum_extensions) + i * n_orders + (order - 1)

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
                        rainbow.fluxlike[c] = this_quantity * 1

                # loop through all the columns in the data extension
                for column in hdu[e].columns:

                    # get a lower case name for the unit
                    c = column.name.lower()

                    # store in the appropriate column of fluxlike array
                    rainbow.fluxlike[c][:, integration_counter - 1] = (
                        hdu[e].data[c] * column_units[c] * 1
                    )

                # increment the running integration total
                integration_counter += 1

    # count up how many times were filled up (and make sure they match)
    n_filled_times = np.sum(np.any(np.isfinite(rainbow.flux), rainbow.waveaxis))
    if n_filled_times != rainbow.ntime:
        cheerfully_suggest(
            f"""
        The x1dints header(s) indicate there should be {rainbow.ntime} integrations,
        but only {n_filled_times} columns of the flux array were populated. Are you
        perhaps missing some segment files? These would be multiple files with
        different numbers after the string `seg` in the filename.
        """
        )

    # try to figure out the errors from available column keywords
    try:
        rainbow.fluxlike["uncertainty"] = rainbow.fluxlike["flux_error"] * 1
    except KeyError:
        rainbow.fluxlike["uncertainty"] = rainbow.fluxlike["error"] * 1

    if "uncertainty" not in rainbow.fluxlike:
        cheerfully_suggest(
            f"""
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
        )
