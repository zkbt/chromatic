"""
Define a reader for NIRISS/SOSS extract1d files produced using the ATOCA
algorithm.
"""
from ...imports import *

__all__ = ["from_atoca"]


def get_time_axis(header):
    t_start = header['DATE-OBS']+'T'+header['TIME-OBS']
    tframe = header['TFRAME']/3600/24
    nframe = header['NFRAMES']
    ngroup = header['NGROUPS']+1
    nint = header['NINTS']
    t_start = Time(t_start, format='isot', scale='utc')
    t = np.arange(nint)*tframe*nframe*ngroup + t_start.jd
    return t*u.d, nint


def from_atoca(rainbow, filepath, order=1):
    if order not in [1, 2]:
        msg = 'Only orders 1 and 2 are extracted by ATOCA.'
        raise NotImplementedError(msg)
    else:
        print('Unpacking order {}'.format(order), flush=True)
    filenames = expand_filenames(filepath)

    for i_file, f in enumerate(tqdm(filenames)):
        hdu_list = fits.open(f)

        if i_file == 0:
            quantities = {}
            times, nints = get_time_axis(hdu_list['PRIMARY'].header)
            rainbow.timelike['time'] = times
            first_time = True

        for i in range(1, len(hdu_list)):
            if hdu_list[i].header['EXTNAME'] != 'EXTRACT1D':
                continue
            if hdu_list[i].header['SPORDER'] != order:
                continue

            for quantity in hdu_list[i].data.names:
                if first_time:
                    quantities[quantity] = hdu_list[i].data[quantity]
                else:
                    quantities[quantity] = np.vstack([quantities[quantity],
                                                      hdu_list[i].data[quantity]])
            first_time = False 

    for quantity in quantities.keys():
        if quantity == 'FLUX':
            rainbow.fluxlike['flux'] = quantities[quantity].T
        elif quantity == 'FLUX_ERROR':
            rainbow.fluxlike['uncertainty'] = quantities[quantity].T
        elif quantity == 'WAVELENGTH':
            rainbow.wavelike['wavelength'] = np.nanmedian(quantities[quantity],
                                                          axis=0)*u.micron
        else:
            rainbow.fluxlike[quantity] = quantities[quantity].T

    n_filled_times = rainbow.fluxlike['flux'].shape[1]
    if n_filled_times != nints:
        warnings.warn(f"""The x1dints header(s) indicate there should be 
            {rainbow.ntime} integrations, but only {n_filled_times} columns of
            the flux array were populated. Are you perhaps missing some 
            segment files?""")
