"""
Define a reader for NIRISS/SOSS extract1d files produced using the ATOCA
algorithm.
"""
from ...imports import *

__all__ = ["from_atoca"]


def get_time_axis(header):
    """Define the time axis for a NIRISS/SOSS TSO

    Parameters
    ----------
    header : Header
        Header of the primary HDU.

    Returns
    -------
    t : array-like
        TSO time axis in jd.
    nint : int
        Number of integrations in the TSO.
    """

    # Each integration is not time stamped - need to reconstruct the time axis
    # based on available information.
    # Get the observation start time
    t_start = header["DATE-OBS"] + "T" + header["TIME-OBS"]
    # Get frame time and convert to days.
    tframe = header["TFRAME"] * u.s.to(u.d)
    # Get number of frames, groups and integrations.
    # nframe willl generally be 1 for NIRISS/SOSS, but just in case...
    nframe = header["NFRAMES"]
    ngroup = header["NGROUPS"] + 1  # Need to add one for reset.
    nint = header["NINTS"]
    # Construct time axis and give units of days.
    t_start = Time(t_start, format="isot", scale="utc")
    t = np.arange(nint) * tframe * nframe * ngroup + t_start.jd
    t *= u.d

    return t, nint


def from_atoca(rainbow, filepath, order=1):
    """Initialize a Rainbow based on the outputs of the ATOCA 1D extraction
    algorithm. This will only be used for NIRISS/SOSS spectra.

    Parameters
    ----------
    rainbow : Rainbow
        The object to be populated. (This is intended to be used as a class
        method of Rainbow or of a class derived from Rainbow, as a way of
        initializing an object from files.)
    filepath : str
        The path to the file, or a glob-friendly string pointing to a group of
        files that should all be loaded together (for example, if an exposure
        was split into multiple segments).
    order : int
        Diffraction order to be unpacked, either 1 or 2. Currently both cannot
        be unpacked simultaneously.
    """

    # Verify that the requested order makes sense.
    if order not in [1, 2]:
        msg = "Only orders 1 and 2 are extracted by ATOCA."
        raise NotImplementedError(msg)
    else:
        # Ensure that the user knows which order they are getting.
        other_order = {1: 2, 2: 1}[order]
        cheerfully_suggest(
            f"""
        You are loading NIRISS order {order} from the ATOCA file
        {filepath}
        If you want the other order, try providing the
        `order={other_order}` as a keyword argument to `Rainbow()`.
        """
        )

    filenames = expand_filenames(filepath)

    # Loop over all input filenames.
    for i_file, f in enumerate(tqdm(filenames, leave=False)):
        hdu_list = fits.open(f)

        # Useful initializations for later.
        quantities = {}
        first_time = True

        # For the first file, define common time axis.
        if i_file == 0:
            # Create time axis
            times, nints = get_time_axis(hdu_list["PRIMARY"].header)
            rainbow.timelike["time"] = times * 1

        # Loop over all extensions.
        for i in range(1, len(hdu_list)):
            # Only consider extract1d extensions of the correct order.
            if hdu_list[i].header["EXTNAME"] != "EXTRACT1D":
                continue
            if hdu_list[i].header["SPORDER"] != order:
                continue
            # Unpack each of the data types stored in each extension.
            # This includes the usual wavelength, flux, and error, but also
            # many other things which may or may not be useful (or even
            # populated). But keep them just in case.
            for quantity in hdu_list[i].data.names:
                if first_time:
                    quantities[quantity] = hdu_list[i].data[quantity]
                else:
                    quantities[quantity] = np.vstack(
                        [quantities[quantity], hdu_list[i].data[quantity]]
                    )
            first_time = False

    # Pack all the above data into a Rainbow object.
    for quantity in quantities.keys():
        # Try to keep the main key names consistant with other formats.
        # Flux-like things in electron/s and wavelengths in microns
        if quantity == "FLUX":
            rainbow.fluxlike["flux"] = quantities[quantity].T * 1.6 * u.electron / u.s
        elif quantity == "FLUX_ERROR":
            rainbow.fluxlike["uncertainty"] = (
                quantities[quantity].T * 1.6 * u.electron / u.s
            )
        elif quantity == "WAVELENGTH":
            rainbow.wavelike["wavelength"] = (
                np.nanmedian(quantities[quantity], axis=0) * u.micron * 1
            )
        else:
            rainbow.fluxlike[quantity] = quantities[quantity].T * 1

    # Warn user if the number of unpacked integrations doesn't match the
    # expected amount.
    n_filled_times = rainbow.fluxlike["flux"].shape[1]
    if n_filled_times != nints:
        cheerfully_suggest(
            f"""The extract1d header(s) indicate there should be
            {rainbow.ntime} integrations, but only {n_filled_times} columns of
            the flux array were populated. Are you perhaps missing some
            segment files?"""
        )
