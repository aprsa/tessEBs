import numpy as np
from astroquery.mast import Catalogs as cat, Observations as obs, Tesscut as tc
from astroquery.simbad import Simbad as simbad


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

    other_ids = simbad.query_objectids(f'TIC {tess_id}')
    if other_ids is not None:
        for other_id in other_ids['ID']:
            asas_id = other_id.split(' ')[1] if 'ASAS' in other_id else None
            kepler_id = other_id.split(' ')[1] if 'Kepler' in other_id else None
            gaia_dr3_id = other_id.split(' ')[2] if 'DR3' in other_id else None
    else:
        asas_id = None
        kepler_id = None
        gaia_dr3_id = None

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
