"""
Define a reader for STScI pipeline x1dints.fits files.
"""
from ...imports import *

__all__ = ["from_x1dints"]


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
    for i_file, f in enumerate(tqdm(filenames)):

        # open this fits file
        hdu = fits.open(f)

        # if this is the first one, populate the shared stuff
        if i_file == 0:

            # grab the entire primary header
            rainbow.metadata["header"] = hdu["PRIMARY"].header

        # set the index to the start of this segment
        integration_counter = hdu["PRIMARY"].header["INTSTART"]
        n_integrations_total = hdu["PRIMARY"].header["NINTS"]
        n_integrations_in_this_segment = (
            hdu["PRIMARY"].header["INTEND"] - hdu["PRIMARY"].header["INTSTART"] + 1
        )

        if hdu["PRIMARY"].header["INSTRUME"] == "NIRISS":
            n_orders = 3
        else:
            n_orders = 1
        assert (order >= 1) and (order <= n_orders)

        non_spectrum_extensions = [
            x for x in ["PRIMARY", "SCI", "INT_TIMES", "ASDF"] if x in hdu
        ]
        if "int_times" in hdu:
            warnings.warn(
                f"""
            It looks like the `x1dints` file you're trying to load
            contains an `int_times` extension, suggesting it's an
            early simulated dataset. If it doesn't succeed with this
            updated reader, please try again with format='x1dints_kludge',
            which was tuned to a bunch of the quirks of ERS Data
            Challenge simulations.
            """
            )

        # make sure sizes match, ignoring PRIMARY, SCI, ASDF
        assert n_integrations_in_this_segment * n_orders == (
            len(hdu) - len(non_spectrum_extensions)
        )

        # loop through the integrations in this segment
        for i in range(n_integrations_in_this_segment):
            e = 2 + i * n_orders + (order - 1)

            # do this stuff only on the very first time through
            if i_file == 0 and i == 0:

                # pull out a wavelength grid (assuming it'll stay constant)
                unit_string = hdu[e].columns["wavelength"].unit
                if unit_string is None:
                    unit_string = "micron"
                    warnings.warn("No wavelength unit was found; assuming micron.")
                wavelength_unit = u.Unit(unit_string)
                rainbow.wavelike["wavelength"] = (
                    hdu[e].data["wavelength"] * wavelength_unit * 1
                )

                # make up a time grid (this is still a little kludgy; it'd be better to have int_times)
                first_mjd_barycentric_integration_midpoint = (
                    hdu["sci"].header["TDB-BEG"] * u.day
                    + hdu["primary"].header["EFFINTTM"] / 2 * u.second
                )
                last_mjd_barycentric_integration_midpoint = (
                    hdu["sci"].header["TDB-END"] * u.day
                    - hdu["primary"].header["EFFINTTM"] / 2 * u.second
                )
                mjd_barycentric_integration_midpoints = np.linspace(
                    first_mjd_barycentric_integration_midpoint,
                    last_mjd_barycentric_integration_midpoint,
                    n_integrations_total,
                )
                astropy_times = Time(
                    mjd_barycentric_integration_midpoints, format="mjd", scale="tdb"
                )
                rainbow.set_times_from_astropy(astropy_times, is_barycentric=True)
                warnings.warn(
                    f"""
                Times were set by linearly interpolating between the exposure
                start and end points. It's very possible these times are off
                by at least a few seconds and possibly up to the duration
                of one integration (= {hdu["primary"].header["EFFINTTM"]}s)
                """
                )

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

        # when done with this file, make sure indices line up
        # assert integration_counter == hdu["PRIMARY"].header["INTEND"]
        #  FIXME - it looks like different instruments may have different INTENDs?
    n_filled_times = np.sum(np.any(np.isfinite(rainbow.flux), rainbow.waveaxis))
    if n_filled_times != rainbow.ntime:
        warnings.warn(
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
        warnings.warn(message)
