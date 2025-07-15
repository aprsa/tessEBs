import os
import numpy as np
import time
import matplotlib.pyplot as plt

from astroquery.mast import Catalogs as cat, Observations as obs, Tesscut as tc
from astroquery.simbad import Simbad as simbad
from astropy.io import fits
from astropy.table import Table

from . import provenances


def bjd2phase(times, bjd0, period, pshift=0.0):
    """
    Convert BJD times to phase.

    Parameters
    ----------
    times : ndarray
        BJD times.
    bjd0 : float
        BJD of the first primary eclipse.
    period : float
        Period of the binary star.
    pshift : float, optional, default=0.0
        Phase shift.

    Returns
    -------
    phases : ndarray
        Phases of the binary star.
    """

    phases = -0.5+((times-bjd0-(pshift+0.5)*period) % period) / period
    return phases


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
        # it seems that id column can be either 'id' or 'ID',
        # depending on the version of astroquery. Let's handle that:
        if 'id' in other_ids.colnames:
            other_ids = other_ids['id']
        elif 'ID' in other_ids.colnames:
            other_ids = other_ids['ID']
        else:
            other_ids = []

        for other_id in other_ids:
            if 'ASAS' in other_id:
                asas_id = other_id.split(' ')[1]
            if 'KIC' in other_id:
                kepler_id = other_id.split(' ')[1]
            if 'DR3' in other_id:
                gaia_dr3_id = other_id.split(' ')[2]

    sectors = [int(s) for s in tc.get_sectors(objectname=f'TIC {tess_id}')['sector'].data]
    provenances = [
        str(provenance) for provenance in set(obs.query_criteria(target_name=tess_id, dataproduct_type='timeseries', project='TESS')['provenance_name'])
    ]

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


def load_data(tess_id, datatype='lc', provenance=None, data_dir='static/catalog/lc_data'):
    """
    Load data from a FITS file.

    Parameters
    ----------
    tess_id : int
        TESS ID of the target.
    datatype : str
        Data type of the target (e.g. 'lc', 'spd').
    provenance : str
        Provenance name. If None, all provenances are used.
    data_dir : str
        Directory where the FITS file is located.

    Returns
    -------
    data : dict
        Dictionary containing data arrays with provenance names as keys.

    Raises
    ------
    FileNotFoundError
        If the FITS file does not exist.
    ValueError
        If the provenance is not found in the FITS file.
        If the data type is not recognized.
    """

    filename = f'{data_dir}/tic{tess_id:010d}.fits'
    if not os.path.exists(filename):
        raise FileNotFoundError(f'File {filename} not found.')

    result = {}
    with fits.open(f'{data_dir}/tic{tess_id:010d}.fits') as hdul:
        # get all provenances in the fits file:
        provenances = list(hdul[0].header['PROVS'].split(','))

        # if provenance is given, use it; otherwise include all provenances:
        if provenance is not None:
            if provenance not in provenances:
                raise ValueError(f'Provenance {provenance} not found in file {filename}.')

            if datatype == 'lc':
                result[provenance] = hdul[provenance].data
            elif datatype == 'spd':
                result[provenance] = hdul[provenance+'-SPD'].data
            else:
                raise ValueError(f'load_data(): data type {datatype} not recognized.')
        else:
            for prov in provenances:
                if datatype == 'lc':
                    result[prov] = hdul[prov].data
                elif datatype == 'spd':
                    result[prov] = hdul[prov+'-SPD'].data
                else:
                    raise ValueError(f'load_data(): data type {datatype} not recognized.')

        return result


def plot_plc(filename, data, bjd0, period, pshift=0.0, title=None, xlabel='Phase', ylabel='Normalized Flux'):
    """
    Plot the phased light curve.

    Parameters
    ----------
    filename : str
        Filename to save the plot.
    data : ndarray
        Data array.
    bjd0 : float
        BJD of the first primary eclipse.
    period : float
        Period of the binary star.
    pshift : float, optional, default=0.0
        Phase shift.
    title : str, optional, default=None
        Title of the plot.
    xlabel : str, optional, default='Phase'
        X-axis label.
    ylabel : str, optional, default='Normalized Flux'
        Y-axis label.

    Returns
    -------
    None
    """

    plt.figure('plc', figsize=(8, 5))
    plt.xlabel(xlabel or 'Phase')
    plt.ylabel(ylabel or 'Flux')
    plt.title(title or '')

    # Get a color cycle from matplotlib
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    for i, (provenance, prov_data) in enumerate(data.items()):
        color = colors[i % len(colors)]
        times, fluxes = prov_data['times'], prov_data['fluxes']
        phases = bjd2phase(times, bjd0, period, pshift=pshift)

        # Plot all points with the same color for consistency
        plt.plot(phases, fluxes, '.', label=provenance, color=color)
        plt.plot(phases[phases > 0.4]-1.0, fluxes[phases > 0.4], '.', color=color)
        plt.plot(phases[phases < -0.4]+1.0, fluxes[phases < -0.4], '.', color=color)

    plt.legend()
    plt.savefig(filename)
    plt.close()
