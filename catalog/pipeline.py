import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import tempfile
import os
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares as BLS
from astroquery.mast import Observations as O

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


def download_fits(tess_id, local_path=None):
    all_data = O.query_object(f'TIC {tess_id}')
    tdp = all_data[(all_data['obs_collection'] == 'TESS') & (all_data['dataproduct_type'] == 'timeseries') & (all_data['target_name'] == f'{tess_id}')]
    for entry in tdp:
        if 'dvt' in entry['dataURL']:
            continue
        print('downloading to ' + os.path.join(local_path, os.path.basename(entry['dataURL'])))
        O.download_file(entry['dataURL'], local_path=os.path.join(local_path, os.path.basename(entry['dataURL'])))


def read_from_all_fits(tess_id, normalize=True, filelist=None):
    if filelist is None:
        fits_list = glob.glob('%s/**/*%d*lc.fits' % (DATADIR, tess_id), recursive=True)
    else:
        fits_list = [x for x in filelist if '%d' % tess_id in x]

    if len(fits_list) == 0:
        return None, None, None

    times, fluxes, ferrs = [], [], []
    for filename in fits_list:
        with fits.open(filename) as hdu:
            # pylint: disable=E1101
            time = hdu[1].data['TIME']
            flux = hdu[1].data['PDCSAP_FLUX']
            flux_err = hdu[1].data['PDCSAP_FLUX_ERR']

            cond = np.where(np.logical_and(~np.isnan(time), ~np.isnan(flux), ~np.isnan(flux_err)))
            time = time[cond]
            if normalize:
                flux_err = flux_err[cond]/flux[cond].mean()
                flux = flux[cond]/flux[cond].mean()
            else:
                flux_err = flux_err[cond]
                flux = flux[cond]
        times.append(time)
        fluxes.append(flux)
        ferrs.append(flux_err)

    # Concatenate all arrays:
    times = np.concatenate(times)
    fluxes = np.concatenate(fluxes)
    ferrs = np.concatenate(ferrs)

    # As globbing fits files might be out of order, we need to sort the arrays:
    s = np.argsort(times)

    return (times[s], fluxes[s], ferrs[s])


def read_from_fits(tess_id, normalize=True):
    filename = glob.glob('%s/*%d*fits' % (DATADIR, tess_id))[0]

    with fits.open(filename) as hdu:
        # pylint: disable=E1101
        time = hdu[1].data['TIME']
        flux = hdu[1].data['PDCSAP_FLUX']
        flux_err = hdu[1].data['PDCSAP_FLUX_ERR']

        cond = np.where(np.logical_and(~np.isnan(time), ~np.isnan(flux), ~np.isnan(flux_err)))
        time = time[cond]
        if normalize:
            flux = flux[cond]/flux[cond].mean()
            flux_err = flux_err[cond]/flux.mean()
        else:
            flux = flux[cond]
            flux_err = flux_err[cond]

    return (time, flux, flux_err)


def bjd2phase(times, bjd0, period, pshift=0.0):
    phases = -0.5+((times-bjd0-(pshift+0.5)*period) % period) / period
    return phases


def run_lombscargle(time, flux, ferr, pmin=0.1, pmax=15, subsample=0.001, npeaks=3, extras=''):
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
