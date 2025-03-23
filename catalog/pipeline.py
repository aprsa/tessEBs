import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os

from astropy.io import fits
from astropy.timeseries import BoxLeastSquares as BLS

from scipy.signal import lombscargle as LS
from scipy.signal import find_peaks, peak_prominences
from scipy.stats import exponweib


def fits_exists(tess_id):
    fits_list = glob.glob('%s/**/*%d*lc.fits' % (DATADIR, tess_id), recursive=True)
    if len(fits_list) > 0:
        return True
    else:
        return False


def read_from_path(tess_id, provenance, data_dir='static/catalog', normalize=True, remove_nans=True):
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
    # logfap = np.log10(1.0 - exponweib.cdf(power, a=1.0, c=1.0, loc=0.0, scale=1.0))

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
