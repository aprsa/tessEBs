import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os

from astropy.io import fits
from astropy.timeseries import BoxLeastSquares as BLS
from astroquery.mast import Observations as O

from scipy.signal import lombscargle as LS
from scipy.signal import find_peaks, peak_prominences
from scipy.stats import exponweib

DATADIR = '/home/users/andrej/projects/tess/data'


def scan_directory(datadir, sector_cap=None):
    fits_list = np.array(glob.glob('%s/**/*lc.fits' % (datadir), recursive=True), dtype=str)
    if sector_cap is not None:
        sectors = np.array([int(f.split('-')[1][1:]) for f in fits_list], dtype=int)
        fits_list = fits_list[sectors <= sector_cap]
    if len(fits_list) == 0:
        print('No fits files found in the %s directory.' % (datadir))
    return fits_list


def fits_exists(tess_id):
    fits_list = glob.glob('%s/**/*%d*lc.fits' % (DATADIR, tess_id), recursive=True)
    if len(fits_list) > 0:
        return True
    else:
        return False


def download_fits(tess_id, destination=None, **kwargs):
    """
    Download fits files for a given TESS ID.

    Keyword arguments:
    ------------------
    obs_collection: str
        Observation collection, 'TESS' or 'HLSP'. By default all collections are downloaded.
    provenance_name: str
        Provenance name. By default all provenance names are downloaded.
    """

    data = O.query_criteria(target_name=tess_id, dataproduct_type='timeseries', project='TESS', **kwargs)
    if len(data) > 0:
        return O.download_products(O.get_product_list(data), download_dir=destination)
    return None


def read_from_path(tess_id, provenance, data_dir=DATADIR, normalize=True, remove_nans=True):
    times, fluxes, ferrs = provenance.lc(tess_id=tess_id, data_dir=data_dir, normalize=normalize, remove_nans=remove_nans)
    return times, fluxes, ferrs


def load_data(tess_id, datatype='lc', provenance=None, data_dir='static/catalog/lc_data'):
    filename = f'{data_dir}/tic{tess_id:010d}.fits'
    if not os.path.exists(filename):
        raise FileNotFoundError(f'File {filename} not found.')

    with fits.open(f'{data_dir}/tic{tess_id:010d}.fits') as hdul:
        if provenance is None:
            provenance = hdul[1].header['EXTNAME']

        if datatype == 'lc':
            return hdul[provenance].data

        elif datatype == 'spd':
            return hdul[provenance+'-SPD'].data

        else:
            raise ValueError(f'Data type {datatype} not recognized.')


def bjd2phase(times, bjd0, period, pshift=0.0):
    phases = -0.5+((times-bjd0-(pshift+0.5)*period) % period) / period
    return phases


def estimate_period(time, flux, window):
    """
    Estimate period of a light curve using the autocorrelation function.

    Parameters
    ----------
    time : ndarray
        Time array.
    flux : ndarray
        Flux array.
    window : float
        Window size for the template.

    Returns
    -------
    period : float
        Estimated period.
    """
    dtime = np.diff(time)
    gaps = np.where(dtime > 1.1 * np.min(dtime))[0]  # not used at present, but it should be used to remove gaps
    dflux = np.diff(flux)
    segment = time[:-1] < time[0] + window
    template = dflux[segment]
    acorr = np.correlate(dflux, template, mode='valid')
    peaks = find_peaks(acorr, height=0.0)
    prominences, _, _ = peak_prominences(acorr, peaks[0])
    peaks = find_peaks(acorr, height=0.0, prominence=np.mean(prominences))
    period = 2 * np.mean(np.diff(peaks[0])) * dtime[0]
    return period


def run_lombscargle(data, pmin=0.1, pmax=15, pstep=0.001, npeaks=0):
    wmin = 2*np.pi/pmax
    wmax = 2*np.pi/pmin
    wstep = 2*np.pi*pstep/pmin/pmax

    w = np.arange(wmin, wmax, wstep)

    time = data['times']
    flux = (data['fluxes'] - np.mean(data['fluxes']))/np.std(data['fluxes'])
    power = LS(time, flux, w, normalize=False)
    logfap = np.log10(1.0 - exponweib.cdf(power, a=1.0, c=1.0, loc=0.0, scale=1.0))

    # Find peaks
    if npeaks > 0:
        peak_indices = (-power).argsort()[:npeaks]
        peak_freqs = w[peak_indices]
        peak_powers = power[peak_indices]
        peak_periods = 1.0/peak_freqs
    else:
        peak_freqs = []
        peak_powers = []
        peak_periods = []

    return {
        'freqs': w,
        'periods': 2*np.pi/w,
        'powers': power,
        'logfaps': logfap,
        'peak_frequencies': peak_freqs,
        'peak_powers': peak_powers,
        'peak_periods': peak_periods
    }


def run_vartools_lombscargle(time, flux, ferr, pmin=0.1, pmax=15, subsample=0.001, npeaks=3, extras=''):
    with tempfile.NamedTemporaryFile() as lcf:
        np.savetxt(lcf, np.vstack((time, flux, ferr)).T)

        cmdline = 'vartools -i %s -ascii -LS %f %f %f %d 1 /%s/ %s' % (lcf.name, pmin, pmax, subsample, npeaks, tempfile.gettempprefix(), extras)
        ls = subprocess.run(cmdline.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        lsres = ls.stdout.split()

        freq, ls, logfap = np.loadtxt('%s.ls' % lcf.name, unpack=True)
        os.unlink('%s.ls' % lcf.name)

        rv = {}
        rv['freq'] = freq
        rv['ls'] = ls
        rv['logfap'] = logfap
        rv['periods'] = [float(lsres[1+4*i+0]) for i in range(npeaks)]
        rv['logfaps'] = [float(lsres[1+4*i+1]) for i in range(npeaks)]
        rv['persnrs'] = [float(lsres[1+4*i+2]) for i in range(npeaks)]
        rv['lsstats'] = [float(lsres[1+4*i+3]) for i in range(npeaks)]

    return rv


def run_bls(time, flux, ferr, pmin=0.1, duration=0.05):
    bls = BLS(t=time, y=flux, dy=ferr)
    periodogram = bls.autopower(duration, minimum_n_transit=2, minimum_period=pmin, maximum_period=(time[-1]-time[0])/2)
    max_power = np.argmax(periodogram.power)
    stats = bls.compute_stats(periodogram.period[max_power], periodogram.duration[max_power], periodogram.transit_time[max_power])
    return periodogram, stats


def process_lc(tess_id, signal_id, bjd0=None, period=None):
    time, flux, flux_err = read_from_fits(tess_id)

    # light curve:
    plt.title('TIC %010d.%02d   t0=%5.5f   P0=%5.5f' % (tess_id, signal_id, bjd0, period))
    plt.xlabel('Time')
    plt.ylabel('Normalized flux')
    plt.plot(time, flux, 'b.')
    plt.savefig('lcs/tic%010d.%02d.lc.png' % (tess_id, signal_id))
    plt.clf()

    # run period finder(s):
    # lsres = run_lombscargle(time, flux, flux_err, subsample=0.001)

    # plot periodogram(s):
    # plt.xlabel('Period')
    # plt.ylabel('LS')
    # plt.plot(1./lsres['freq'], lsres['ls'], 'b-')
    # for i in range(3):
    #     plt.axvline(lsres['periods'][i], ls='--', c='r')
    # if period is not None:
    #     plt.axvline(period, ls='-', c='g')
    # plt.show()

    phase = bjd2phase(time, bjd0=bjd0, period=period)

    plt.title('TIC %010d.%02d   t0=%5.5f   P0=%5.5f' % (tess_id, signal_id, bjd0, period))
    plt.xlabel('Phase')
    plt.ylabel('Normalized flux')
    plt.plot(phase, flux, 'b.')
    plt.savefig('lcs/tic%010d.%02d.ph.png' % (tess_id, signal_id))
    plt.clf()


def generate_plots(tess_id, prefix, tlc=False, plc=False, spd=False, signal_id=1, bjd0=None, period=None, filelist=None):
    # if ascii verson of the data exists, load it:
    norm_lc = f'catalog/static/catalog/lcs/tic{tess_id:010d}.norm.lc'
    try:
        time, flux, ferr = np.loadtxt(norm_lc, usecols=(0, 1, 2), unpack=True)
    except:
        time, flux, ferr = read_from_all_fits(tess_id, normalize=True, filelist=filelist)

    plt.rcParams.update({'font.family': 'serif', 'font.size': 20})
    plt.figure(figsize=(16,8))

    if tlc:
        plt.xlabel('Time [days]')
        plt.ylabel('Median-normalized flux')
        plt.plot(time, flux, 'b-')
        plt.savefig(f'{prefix}/tic{tess_id:010d}.tlc.png')

    if plc and period is not None:
        plt.clf()
        if bjd0 is None:
            bjd0 = 0.0
        phase = bjd2phase(time, bjd0, period)
        plt.xlabel('Phase')
        plt.ylabel('Median-normalized flux')
        plt.plot(phase, flux, 'b.')
        plt.plot(phase[phase>=0.4]-1, flux[phase>=0.4], 'b.')
        plt.plot(phase[phase<=-0.4]+1, flux[phase<=-0.4], 'b.')
        plt.savefig(f'{prefix}/tic{tess_id:010d}.{signal_id:02d}.plc.png')
    
    if spd:
        try:
            ls = run_lombscargle(time, flux, ferr, pmin=0.1, pmax=1.2*period, npeaks=1)

            plt.clf()
            plt.xlabel('Period [days]')
            plt.ylabel('Spectral power density')
            plt.plot(1/ls['freq'], ls['ls'], 'b-')
            plt.axvline(period, ls='--', color='r')
            plt.savefig(f'{prefix}/tic{tess_id:010d}.{signal_id:02d}.spd.png')
        except:
            print(f'SPD computation for TIC {tess_id:010d} failed; appending to failed.log.')
            with open('failed.log', 'a') as f:
                f.write(f'{tess_id:010d} SPD\n')

    plt.close()

    return
