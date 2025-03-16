import os
import numpy as np
import time
from astroquery.mast import Catalogs as cat, Observations as obs, Tesscut as tc
from astroquery.simbad import Simbad as simbad
from astropy.io import fits
from astropy.table import Table

from . import provenances


def download_meta(tess_id):
    """
    Downloads metadata for a given TESS ID from MAST and SIMBAD.

    Parameters
    ----------
    tess_id : int
        TESS ID of the target.

    Returns
    -------
    dict
        Dictionary with metadata
    """

    mast_entry = cat.query_criteria(catalog='TIC', ID=tess_id)
    if len(mast_entry) == 0:
        return {
            'status': 'error',
            'message': 'No entry found in MAST database.'
        }

    def get(entry, key):
        return entry[key][0] if ~np.isnan(entry[key][0]) else None

    ra = get(mast_entry, 'ra')
    dec = get(mast_entry, 'dec')
    glon = get(mast_entry, 'gallong')
    glat = get(mast_entry, 'gallat')
    Tmag = get(mast_entry, 'Tmag')
    teff = get(mast_entry, 'Teff')
    logg = get(mast_entry, 'logg')
    abun = get(mast_entry, 'MH')
    pmra = get(mast_entry, 'pmRA')
    pmdec = get(mast_entry, 'pmDEC')

    gaia_id = int(mast_entry['GAIA'][0]) if mast_entry['GAIA'].tolist()[0] is not None else None

    # defaults:
    asas_id = None
    kepler_id = None
    gaia_dr3_id = None

    other_ids = simbad.query_objectids(f'TIC {tess_id}')
    if other_ids is not None:
        for other_id in other_ids['id']:
            if 'ASAS' in other_id:
                asas_id = other_id.split(' ')[1]
            if 'Kepler' in other_id:
                kepler_id = other_id.split(' ')[1]
            if 'DR3' in other_id:
                gaia_dr3_id = other_id.split(' ')[2]

    sectors = [int(s) for s in tc.get_sectors(objectname=f'TIC {tess_id}')['sector'].data]
    provenances = [str(provenance) for provenance in set(obs.query_criteria(target_name=tess_id, dataproduct_type='timeseries', project='TESS')['provenance_name'])]

    meta = {
        'status': 'success',
        'ra': ra,
        'dec': dec,
        'glon': glon,
        'glat': glat,
        'Tmag': Tmag,
        'teff': teff,
        'logg': logg,
        'abun': abun,
        'pmra': pmra,
        'pmdec': pmdec,
        'gaia_id': gaia_id,
        'asas_id': asas_id,
        'kepler_id': kepler_id,
        'gaia_dr3_id': gaia_dr3_id,
        'sectors': sectors,
        'provenances': provenances,
    }

    # convert numpy types to native for json serialization:
    for key, value in meta.items():
        if isinstance(value, np.generic):
            meta[key] = value.item()

    return meta


def download_data(tess_id, dest_dir='static/catalog', **kwargs):
    """
    Download fits files for a given TESS ID.

    Arguments:
    ----------
    tess_id: int, required
        TESS ID of the target.
    dest_dir: str, optional, default='static/catalog'
        Destination directory for the downloaded files.

    Keyword arguments:
    ------------------
    obs_collection: str
        Observation collection, 'TESS' or 'HLSP'. By default all collections are downloaded.
    provenance_name: str
        Provenance name. By default all provenance names are downloaded.

    Returns:
    --------
    list or None
        List of downloaded files. None if no files were found.
    """

    data = obs.query_criteria(target_name=tess_id, dataproduct_type='timeseries', project='TESS', **kwargs)
    if len(data) > 0:
        return obs.download_products(obs.get_product_list(data), download_dir=dest_dir)

    return None


def syndicate_data(tic, dest_dir='static/catalog/lc_data', force_overwrite=False):
    fname = f'{dest_dir}/tic{tic.tess_id:010d}.fits'
    if os.path.exists(fname) and not force_overwrite:
        return fname

    # the method assumes that provenances are up to date.
    if tic.provenances.count() == 0:
        raise ValueError(f'no provenances found for TIC {tic.tess_id:010d}')

    header = fits.Header()
    header['timestmp'] = time.ctime()  # for versioning purposes
    header['provs'] = ','.join([prov.name for prov in tic.provenances.all()])
    hdus = [fits.PrimaryHDU(header=header),]

    for provenance in tic.provenances.all():
        provenance = provenances.get_provenance(provenance.name)
        provenance.download(tess_id=tic.tess_id)
        data = provenance.lc(tic=tic)
        lc = np.array(data, dtype=np.float64).T
        hdus.append(fits.table_to_hdu(Table(lc, names=('times', 'fluxes', 'ferrs'), meta={'extname': provenance.name})))

    fits.HDUList(hdus).writeto(fname, overwrite=True)

    return fname
