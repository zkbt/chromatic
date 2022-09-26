"""
Define a reader for STScI pipeline x1dints.fits files.
"""
from ...imports import *

__all__ = ["from_x1dints_kludge"]


def setup_integration_times(filenames):
    """
    Figure out the times of each integration,
    either by reading them from the headers
    or by totally just making them up from scratch.

    Parameters
    ----------
    filenames : list
        A sorted list of filenames, corresponding to one or more segments.
    """

    # create an empty timelike dictionary
    timelike = {}

    # check the first segment
    first_hdu = fits.open(filenames[0])

    # does the first segment define some times?
    try:
        assert len(first_hdu["int_times"].data) > 0

        # grab the integration time information
        for c in first_hdu["int_times"].data.columns.names:
            timelike[c] = first_hdu["int_times"].data[c]

        # be sure to set our standard time axis
        timelike["time"] = timelike["int_mid_BJD_TDB"] * u.day
    # if times are not in header, make up some imaginary ones!
    except:
        # alert the user to what we're doing
        cheerfully_suggest("No times found! Making up imaginary ones!")
        last_hdu = fits.open(filenames[-1])

        # figure out the total number of integrations (DOES THIS NEED THE -1?)
        try:
            N_integrations = last_hdu["PRIMARY"].header["INTEND"] - 1
        except KeyError:
            # (this kludge necessary for CV-simulated NIRSpec)
            N_integrations = last_hdu["PRIMARY"].header["NINTS"]

        # get the time per integration (DOES THIS INCLUDE OVERHEADS?)
        time_per_integration = last_hdu["PRIMARY"].header["EFFINTTM"] * u.s
        cheerfully_suggest(
            f"The imaginary times assume {time_per_integration}/integration."
        )

        # create a fake array of times
        fake_times = np.arange(N_integrations) * time_per_integration
        timelike["time"] = fake_times.to(u.day)

    # return a timelike dictionary
    return timelike


def from_x1dints_kludge(rainbow, filepath, **kw):
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

    # figure out the times (necessary to set up arrays)
    timelike = setup_integration_times(filenames)
    rainbow.timelike.update(**timelike)

    # loop over file (each one a segment)
    for i_file, f in enumerate(tqdm(filenames, leave=False)):

        # open this fits file
        hdu = fits.open(f)

        # if this is the first one, populate the shared stuff
        if i_file == 0:

            # grab the entire primary header
            rainbow.metadata["header"] = hdu["PRIMARY"].header
            # TO-DO: (watch out for the things that aren't constant across segments!)

        # set the index to the start of this segment
        try:
            integration_counter = hdu["PRIMARY"].header["INTSTART"]
        except KeyError:
            # (this kludge necessary for CV-simulated NIRSpec)
            integration_counter = 0

        try:
            n_integrations_predicted_by_header = (
                hdu["PRIMARY"].header["INTEND"] - hdu["PRIMARY"].header["INTSTART"]
            )
        except KeyError:
            n_integrations_predicted_by_header = 0
        n_integrations_in_segment = np.maximum(
            len(hdu) - 3, n_integrations_predicted_by_header
        )

        # print(integration_counter, n_integrations_in_segment, len(hdu))

        # loop through the spectra
        for e in range(2, 2 + n_integrations_in_segment):

            # do this stuff only on the very first time through
            if i_file == 0 and e == 2:
                # pull out a wavelength grid (assuming it'll stay constant)
                unit_string = hdu[e].columns["wavelength"].unit
                if unit_string is None:
                    unit_string = ""
                wavelength_unit = u.Unit(unit_string)
                rainbow.wavelike["wavelength"] = (
                    hdu[e].data["wavelength"] * wavelength_unit * 1
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
        A long-term solution might be to fix the `from_x1dints_kludge`
        reader, but a short-term solution would be to pick one
        of the columns listed above and say something like

        x.fluxlike['uncertainty'] = x.fluxlike['some-other-relevant-error-column']

        where `x` is the Rainbow you just created.
        """
        cheerfully_suggest(message)

    # try to guess wscale (and then kludge and call it linear)
    # rainbow._guess_wscale()
    # rainbow.metadata['wscale'] = 'linear' # TODO: fix this kludge
