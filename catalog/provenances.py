import glob
import numpy as np

from astroquery.mast import Observations
from astropy.io import fits


class Provenance:
    name = 'generic'
    fits_extension = 1
    object_id = 'tess_id'
    suffix = 'lc.fits'

    def __init__(self):
        pass

    def download(self, tess_id, download_dir='static/catalog'):
        obs = Observations.query_criteria(target_name=tess_id, **self.query_filter)
        if len(obs) == 0:
            return None

        product_list = Observations.get_product_list(obs)

        # remove non-LC products:
        non_lc = [i for i, fn in enumerate(product_list['dataURI']) if self.suffix not in fn]
        product_list.remove_rows(non_lc)
        if len(product_list) == 0:
            return None

        data = Observations.download_products(product_list, download_dir=download_dir, **self.data_filter)
        return data

    def lc(self, tic, data_dir='static/catalog/mastDownload', **kwargs):
        normalize = kwargs.get('normalize', True)
        remove_nans = kwargs.get('remove_nans', True)

        fits_list = glob.glob(f'{data_dir}/**/*{self.name.lower()}*{getattr(tic, self.object_id)}*{self.suffix}', recursive=True)
        if len(fits_list) == 0:
            return None, None, None

        all_times, all_fluxes, all_ferrs = [], [], []

        for fits_file in fits_list:
            with fits.open(fits_file) as hdul:
                try:
                    data = hdul[self.fits_extension].data
                except Exception as e:
                    print(f'Error reading {fits_file}: {e}')
                    continue
                times = data[self.cols['time']]

                # some provenances have multiple flux/ferr keywords based on the version.
                if type(self.cols['flux']) is list:
                    # try all keywords in order:
                    found = False
                    for kw in self.cols['flux']:
                        if kw in data.columns.names:
                            found = True
                            break

                    if not found:
                        raise ValueError(f'Flux keyword not found in {fits_file}.')

                    fluxes = data[kw]
                else:
                    fluxes = data[self.cols['flux']]

                if type(self.cols['ferr']) is list:
                    # try all keywords in order:
                    found = False
                    for kw in self.cols['ferr']:
                        if kw in data.columns.names:
                            found = True
                            break

                    if not found:
                        raise ValueError(f'Flux error keyword not found in {fits_file}.')

                    ferrs = data[kw]
                else:
                    ferrs = data[self.cols['ferr']]

                # flags are either integers or characters:
                flags = data[self.cols['flag']]
                if type(data[self.cols['flag']]) is np.char.chararray:
                    flags = np.array([0 if f == 'G' else -1 for f in flags])

                # remove bad quality flags:
                times = times[flags == 0]
                fluxes = fluxes[flags == 0]
                ferrs = ferrs[flags == 0]

                if remove_nans:
                    mask = ~((times != times) | (fluxes != fluxes) | (ferrs != ferrs))
                    times = times[mask]
                    fluxes = fluxes[mask]
                    ferrs = ferrs[mask]

                if normalize:
                    mean_flux = np.nanmean(fluxes)
                    ferrs /= mean_flux
                    fluxes /= mean_flux

            all_times.append(times)
            all_fluxes.append(fluxes)
            all_ferrs.append(ferrs)

        all_times = np.concatenate(all_times)
        all_fluxes = np.concatenate(all_fluxes)
        all_ferrs = np.concatenate(all_ferrs)

        # sort by time:
        s = np.argsort(all_times)

        return all_times[s], all_fluxes[s], all_ferrs[s]


class SPOCProvenance(Provenance):
    name = 'SPOC'
    fits_extension = 1

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'SPOC',
    }
    data_filter = {
        'productSubGroupDescription': 'LC',
        'description': 'Light curves',
    }
    cols = {
        'time': 'TIME',
        'flux': 'SAP_FLUX',
        'ferr': 'SAP_FLUX_ERR',
        'flag': 'QUALITY'
    }


class TESSSPOCProvenance(Provenance):
    name = 'TESS-SPOC'
    fits_extension = 1

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'TESS-SPOC',
    }
    data_filter = {
        'description': 'FITS',
    }
    cols = {
        'time': 'TIME',
        'flux': 'PDCSAP_FLUX',
        'ferr': 'PDCSAP_FLUX_ERR',
        'flag': 'QUALITY'
    }


class QLPProvenance(Provenance):
    name = 'QLP'
    fits_extension = 1

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'QLP',
    }
    data_filter = {
        'description': 'FITS',
    }
    cols = {
        'time': 'TIME',
        'flux': ['DET_FLUX', 'KSPSAP_FLUX'],
        'ferr': ['DET_FLUX_ERR', 'KSPSAP_FLUX_ERR'],
        'flag': 'QUALITY'
    }


class GSFCEleanorLiteProvenance(Provenance):
    name = 'GSFC-ELEANOR-LITE'
    fits_extension = 1

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'GSFC-ELEANOR-LITE',
    }
    data_filter = {
    }
    cols = {
        'time': 'TIME',
        'flux': 'CORR_FLUX',
        'ferr': 'FLUX_ERR',
        'flag': 'QUALITY'
    }


class TGLCProvenance(Provenance):
    name = 'TGLC'
    fits_extension = 1
    object_id = 'gaia_id'

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'TGLC',
    }
    data_filter = {
    }
    cols = {
        'time': 'time',
        'flux': 'aperture_flux',
        'ferr': 'background',
        'flag': 'TESS_flags'
    }


class CDIPSProvenance(Provenance):
    name = 'CDIPS'
    fits_extension = 1
    object_id = 'gaia_id'

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'CDIPS',
    }
    data_filter = {
    }
    cols = {
        'time': 'TMID_BJD',
        'flux': 'IFL1',
        'ferr': 'IFE1',
        'flag': 'IRQ1'
    }


class TASOCProvenance(Provenance):
    name = 'TASOC'
    fits_extension = 1

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'TASOC',
    }
    data_filter = {
    }
    cols = {
        'time': 'TIME',
        'flux': 'FLUX_CORR',
        'ferr': 'FLUX_CORR_ERR',
        'flag': 'QUALITY'
    }


class DiamanteProvenance(Provenance):
    name = 'DIAMANTE'
    fits_extension = 1

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'DIAMANTE',
    }
    data_filter = {
    }
    cols = {
        'time': 'BTJD',
        'flux': 'LC1_AP1',
        'ferr': 'ELC1_AP1',
        'flag': 'FLAG_AP1'
    }


class PathosProvenance(Provenance):
    name = 'PATHOS'
    fits_extension = 1

    query_filter = {
        'dataproduct_type': 'timeseries',
        'project': 'TESS',
        'provenance_name': 'PATHOS',
    }
    data_filter = {
    }
    cols = {
        'time': 'TIME',
        'flux': ['PSF_FLUX_COR', 'AP1_FLUX_COR'],
        'ferr': 'SKY_LOCAL',
        'flag': 'DQUALITY'
    }


def get_provenance(name):
    all_provenances = Provenance.__subclasses__()
    for provenance in all_provenances:
        if provenance.name == name:
            return provenance()
    raise ValueError(f'Provenance subclass for {name} not found.')
