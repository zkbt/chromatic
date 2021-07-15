from ...imports import *

__all__ = ['from_x1dints', 'expand_filenames']

def expand_filenames(filepath):
    if type(filepath) == list:
        filenames = filepath
    else:
        filenames = np.sort(glob.glob(filepath))
    return filenames

def from_x1dints(rainbow, filepath):
    '''
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
    '''

    # create a list of filenames (which might be just 1)
    filenames = expand_filenames(filepath)


    # loop over file (each one a segment)
    for i_file, f in enumerate(tqdm(filenames)):

        # open this fits file
        hdu = fits.open(f)

        # if this is the first one, populate the shared stuff
        if i_file == 0:

            # grab the entire primary header
            rainbow.metadata['header'] = hdu['PRIMARY'].header
            # TO-DO: (watch out for the things that aren't constant across segments!)

            # grab the integration time information
            for c in hdu['int_times'].data.columns.names:
                rainbow.timelike[c] = hdu['int_times'].data[c]

            # be sure to set our standard time axis
            rainbow.timelike['time'] = rainbow.timelike['int_mid_BJD_TDB']*u.day

        # set the index to the start of this segment
        i_integration = hdu['PRIMARY'].header['INTSTART']-1

        # loop through the spectra
        for e in range(2,len(hdu)-1):

            # do this stuff only on the very first time through
            if i_file == 0 and e == 2:
                # pull out a wavelength grid (assuming it'll stay constant)
                wavelength_unit = u.Unit(hdu[e].columns['wavelength'].unit)
                rainbow.wavelike['wavelength'] = hdu[e].data['wavelength']*wavelength_unit

                # set up the fluxlike quantities
                column_units = {}
                for column in hdu[e].columns:

                    # get a lower case name for the unit
                    c = column.name.lower()

                    # set up an empty Quantity array, with units if known
                    this_quantity = u.Quantity(np.nan*np.ones((rainbow.nwave, rainbow.ntime)))
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
                rainbow.fluxlike[c][:, i_integration] = hdu[e].data[c]*column_units[c]

            # increment the running integration total
            i_integration += 1

        # when done with this file, make sure indices line up
        assert i_integration == hdu['PRIMARY'].header['INTEND']

    n_filled_times = np.sum(np.any(np.isfinite(rainbow.flux), rainbow.waveaxis))
    if n_filled_times != rainbow.ntime:
        warnings.warn(f'''
        The x1dints header(s) indicate there should be {rainbow.ntime} integrations,
        but only {n_filled_times} columns of the flux array were populated. Are you
        perhaps missing some segment files?
        ''')
    # try to guess wscale (and then kludge and call it linear)
    #rainbow._guess_wscale()
    #rainbow.metadata['wscale'] = 'linear' # TODO: fix this kludge

    # remove units
    rainbow.fluxlike['flux'] /= u.Jy # TODO: fix and/or test this kludge
