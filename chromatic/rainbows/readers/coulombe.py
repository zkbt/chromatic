"""
Define a reader for STScI pipeline x1dints.fits files.
"""
from ...imports import *

__all__ = ["from_coulombe"]


def from_coulombe(rainbow, filepath, order=1):
    """
    Populate a Rainbow from an STScI pipeline x1dints file,
    or a group of x1dints files for multiple segments.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated.
    filepath : str
        The path to the file, or a glob-friendly string
        pointing to a group of files that should all be
        loaded together (for example, if an exposure was
        split into multiple segments).
    order : int
        The spectral order to load.
    """

    # create a list of filenames (which might be just 1)
    filenames = expand_filenames(filepath)

    # figure out the total number of integrations
    last_hdu = fits.open(filenames[-1])
    n_integrations_predicted_by_header = last_hdu["PRIMARY"].header["NINTS"]

    # populate a 1D array of times (with astropy units of time)
    fake_times = np.arange(n_integrations_predicted_by_header) * u.minute
    cheerfully_suggest("The times are totally made up!")
    rainbow.timelike["time"] = fake_times

    # keep track of integrations loaded yet
    integration_counter = 0

    # loop over file (each one a segment)
    for i_file, f in enumerate(tqdm(filenames, leave=False)):

        # open this fits file
        hdu = fits.open(f)

        # there are three orders!
        n_integrations_in_segment = int((len(hdu) - 2) / 3)

        # if this is the first one, populate the shared stuff
        if i_file == 0:
            # grab the entire primary header
            rainbow.metadata["header"] = hdu["PRIMARY"].header

        start_extension = n_integrations_in_segment * (order - 1) + 1
        end_extension = start_extension + n_integrations_in_segment

        # loop through the spectra
        for e in range(start_extension, end_extension):

            # do this stuff only on the very first time through
            if i_file == 0 and e == start_extension:
                # pull out a wavelength grid (assuming it'll stay constant)
                wavelength_unit = u.Unit(hdu[e].columns["wavelength"].unit)
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
                    rainbow.fluxlike[c] = this_quantity * 1

            for column in hdu[e].columns:

                # get a lower case name for the unit
                c = column.name.lower()
                rainbow.fluxlike[c][:, integration_counter - 1] = (
                    hdu[e].data[c] * column_units[c] * 1
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
    other_order = {1: 2, 2: 1}[order]
    message = f"""
    Loading NIRISS order '{order}'. If you want the other order,
    trying `r = Rainbow(..., format='coulombe', order={other_order})`
    """
    cheerfully_suggest(message)

    # try to guess wscale (and then kludge and call it linear)
    # rainbow._guess_wscale()
    # rainbow.metadata['wscale'] = 'linear' # TODO: fix this kludge
