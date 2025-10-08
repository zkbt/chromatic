from astroquery.mast import Observations

__all__ = ["download_from_mast"]


def download_from_mast(
    proposal_id="2734",
    target_name="WASP-96",
    obs_collection="JWST",
    instrument_name="*",
    productSubGroupDescription="X1DINTS",
    calib_level=2,
):
    """
    Download (JWST) data products from MAST.

    Given a set of search criteria, find the data products
    and download them locally. This has been tested only
    minimally for retrieving X1DINTS files for JWST
    time-series observations for quick look analyses.
    It should hopefully work more broadly?

    If the data Exclusive Access data that are still
    in their proprietary period, you will need to
    login with `Observations.login()` to authenticate
    your connection to the MAST archive.

    For more details on programmatically accessing
    data from MAST, it may be useful to explore the
    example notebooks available at:

    https://spacetelescope.github.io/mast_notebooks/

    Parameters
    ----------

    **kw : dict
        Remaining keywords will be passed to `Observations.query_critera`.
        Available keywords for this initial search likely include:
        ['intentType', 'obs_collection', 'provenance_name',
       'instrument_name', 'project', 'filters', 'wavelength_region',
       'target_name', 'target_classification', 'obs_id', 's_ra', 's_dec',
       'proposal_id', 'proposal_pi', 'obs_title', 'dataproduct_type',
       'calib_level', 't_min', 't_max', 't_obs_release', 't_exptime',
       'em_min', 'em_max', 'objID', 's_region', 'jpegURL', 'distance',
       'obsid', 'dataRights', 'mtFlag', 'srcDen', 'dataURL',
       'proposal_type', 'sequence_number']

    Returns
    -------
    downloaded : astropy.table.Table
        A table summarizing what was downloaded. The columns
        of this table will likely include:
        ['Local Path', 'Status', 'Message', 'URL']
        where 'Local Path' indicates to where the file
        was downloaded (relative to the local path),
        and 'Status' indicates whether the download
        was successful.

    """

    # get a table of all observations matching the search criteria
    obs_table = Observations.query_criteria(
        obs_collection=obs_collection,
        proposal_id=proposal_id,
        instrument_name=instrument_name,
        target_name=target_name,
    )

    # get a table of the list of products associated with those observations
    products = Observations.get_product_list(obs_table["obsid"])

    # filter those products down to just what's requested
    if isinstance(calib_level, int):
        calib_level = [calib_level]
    filtered_products = Observations.filter_products(
        products,
        productSubGroupDescription=productSubGroupDescription,
        calib_level=calib_level,
    )

    # download the filtered products
    downloaded = Observations.download_products(filtered_products)

    # return table summarizing the downloaded files
    return downloaded
