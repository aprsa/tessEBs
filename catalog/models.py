from django.db import models
from . import backend
from .pipeline import run_lombscargle, run_bls
from .provenances import get_provenance

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table


class Sector(models.Model):
    sector_id = models.IntegerField('sector id')
    date_start = models.DateField('start date')
    date_end = models.DateField('end date')
    spacecraft_ra = models.FloatField('spacecraft R.A.')
    spacecraft_dec = models.FloatField('spacecraft dec')
    spacecraft_roll = models.FloatField('spacecraft roll')
    camera1_ra = models.FloatField('camera 1 R.A.')
    camera1_dec = models.FloatField('camera 1 dec')
    camera1_roll = models.FloatField('camera 1 roll')
    camera2_ra = models.FloatField('camera 2 R.A.')
    camera2_dec = models.FloatField('camera 2 dec')
    camera2_roll = models.FloatField('camera 2 roll')
    camera3_ra = models.FloatField('camera 3 R.A.')
    camera3_dec = models.FloatField('camera 3 dec')
    camera3_roll = models.FloatField('camera 3 roll')
    camera4_ra = models.FloatField('camera 4 R.A.')
    camera4_dec = models.FloatField('camera 4 dec')
    camera4_roll = models.FloatField('camera 4 roll')

    def __str__(self):
        return '%d' % self.sector_id

    def __repr__(self):
        return (
            f"<Sector(sector_id={self.sector_id}, "
            f"date_start={self.date_start}, "
            f"date_end={self.date_end}, "
            f"spacecraft_ra={self.spacecraft_ra}, "
            f"spacecraft_dec={self.spacecraft_dec}, "
            f"spacecraft_roll={self.spacecraft_roll}, "
            f"camera1_ra={self.camera1_ra}, "
            f"camera1_dec={self.camera1_dec}, "
            f"camera1_roll={self.camera1_roll}, "
            f"camera2_ra={self.camera2_ra}, "
            f"camera2_dec={self.camera2_dec}, "
            f"camera2_roll={self.camera2_roll}, "
            f"camera3_ra={self.camera3_ra}, "
            f"camera3_dec={self.camera3_dec}, "
            f"camera3_roll={self.camera3_roll}, "
            f"camera4_ra={self.camera4_ra}, "
            f"camera4_dec={self.camera4_dec}, "
            f"camera4_roll={self.camera4_roll})>"
        )


class Provenance(models.Model):
    name = models.CharField(max_length=32)

    def __str__(self):
        return f'{self.name}'

    def __repr__(self):
        return (
            f"<Provenance(name='{self.name}')>"
        )


class TIC(models.Model):
    # this class defines the fields from the TESS input catalog:
    tess_id = models.BigIntegerField('tess id')
    ra = models.FloatField('ra', null=True, blank=True)
    dec = models.FloatField('dec', null=True, blank=True)
    glon = models.FloatField('galactic longitude', null=True, blank=True)
    glat = models.FloatField('galactic latitude', null=True, blank=True)
    Tmag = models.FloatField('tess magnitude', null=True, blank=True)
    teff = models.FloatField('effective temperature', null=True, blank=True)
    logg = models.FloatField('surface gravity', null=True, blank=True)
    abun = models.FloatField('metallicity', null=True, blank=True)
    pmra = models.FloatField('proper motion in ra [mas/yr]', null=True, blank=True)
    pmdec = models.FloatField('proper motion in dec [mas/yr]', null=True, blank=True)
    gaia_id = models.BigIntegerField('gaia id', null=True, blank=True)
    gaia_dr3_id = models.BigIntegerField('gaia dr3 id', null=True, blank=True)
    kepler_id = models.BigIntegerField('kepler id', null=True, blank=True)
    asas_id = models.CharField(max_length=20, null=True, blank=True)
    sectors = models.ManyToManyField(Sector, blank=True)
    provenances = models.ManyToManyField(Provenance, blank=True)

    # OBSOLETE:
    class DataType(models.TextChoices):
        FFI = 'FFI', 'Full Frame Image (FFI)'
        TPF = 'TPF', 'Target Pixel File (TPF)'
        LCF = 'LCF', 'Light Curve File (LCF)'
    datatype = models.CharField(max_length=3, choices=DataType.choices, default=DataType.LCF)

    class Meta:
        verbose_name_plural = 'TICs'

    @classmethod
    def from_mast(cls, tess_id, **kwargs):
        """
        Pulls the TIC data from MAST and creates a new instance of the TIC class.
        If the TIC is already in the database, it will update the information
        from MAST.

        Parameters
        ----------
        tess_id : int, required
            TESS ID of the target star.
        syndicate_data : bool, optional, default True
            If True, the method will download the data from MAST.
        attach_spd_data : bool, optional, default True
            If True, the method will attach spectral power density data
            to the syndicated data. Note that `syndicate_data` must also
            be set to True.
        create_static : bool, optional, default False
            If True, the method will create static files for the target.
        static_dir : str, optional, default 'static/catalog'
            The directory where the static files will be stored.
        overwrite_static_files : bool, optional, default False
            If True, the method will overwrite existing static files.

        Returns
        -------
        TIC instance
            An instance of the TIC class with the data from MAST.

        Raises
        ------
        ValueError
            If the TIC is not found in the MAST database.
            If the sector is not found in the internal database.
        """

        meta = backend.download_meta(tess_id)
        status = meta.pop('status')
        if status != 'success':
            raise ValueError(f'TIC {tess_id} not found in the MAST database.')

        sectors = meta.pop('sectors')
        provenances = meta.pop('provenances')

        # if the TIC is already in the database, use it; otherwise instantiate a new one:
        tic, created = cls.objects.get_or_create(tess_id=tess_id, defaults=meta)

        # if new TIC, save it to assign an ID:
        if created:
            tic.save()

        # add sectors:
        for sector_id in sectors:
            sobj = Sector.objects.get(sector_id=sector_id)
            if sobj is None:
                raise ValueError(f'sector {sector_id} not found in the database.')
            tic.sectors.add(Sector.objects.get(sector_id=sector_id))

        # add provenances:
        for pname in provenances:
            prov, _ = Provenance.objects.get_or_create(name=pname)
            tic.provenances.add(prov)

        # download and syndicate data if requested:
        if kwargs.get('syndicate_data', True):
            force_overwrite = kwargs.get('overwrite_static_files', False)
            for provenance_name in provenances:
                provenance = get_provenance(provenance_name)
                provenance.download(tess_id=tess_id)
            tic.syndicate_data(force_overwrite=force_overwrite)

            # attach SPD data if requested:
            if kwargs.get('attach_spd_data', True):
                force_overwrite = kwargs.get('overwrite_static_files', True)
                tic.attach_spds_to_data(force_overwrite=force_overwrite)

        # create static files if requested:
        if kwargs.get('create_static', False):
            force_overwrite = kwargs.get('overwrite_static_files', True)
            tic.create_static_files(force_overwrite=force_overwrite)

        return tic

    def syndicate_data(self, data_dir='static/catalog/lc_data', force_overwrite=False):
        return backend.syndicate_data(self, dest_dir=data_dir, force_overwrite=force_overwrite)

    def attach_spds_to_data(self, data_dir='static/catalog/lc_data', force_overwrite=False):
        # the method assumes that provenances are up to date.
        if self.provenances.count() == 0:
            raise ValueError(f'no provenances found for TIC {self.tess_id}')

        with fits.open(f'{data_dir}/tic{self.tess_id:010d}.fits', mode='update') as hdul:
            for provenance in [p.name for p in self.provenances.all()]:
                if provenance+'-SPD' in hdul and not force_overwrite:
                    continue

                spd = run_lombscargle(hdul[provenance].data)
                spd = np.vstack((spd['periods'], spd['powers'])).T

                # remove any old SPD data:
                if provenance+'-SPD' in hdul:
                    hdul.pop(provenance+'-SPD')

                hdul.append(fits.table_to_hdu(Table(spd, names=('periods', 'powers'), meta={'extname': provenance+'-SPD'})))

    def get_lightcurve(self, provenance=None):
        """
        Get the syndicated lightcurve data for the TIC.

        Arguments:
        ---------
        provenance : str
            The provenance name to select the data from. If None, the
            provenance will default to the first one in the list.
        """
        return backend.load_data(self.tess_id, provenance=provenance)

    def get_spd(self, provenance=None):
        """
        Get the syndicated spectral power density data for the TIC.

        Arguments:
        ---------
        provenance : str
            The provenance name to select the data from. If None, the
            provenance will default to the first one in the list.
        """
        return backend.load_data(self.tess_id, datatype='spd', provenance=provenance)

    def run_bls(self, pmin=0.5, pmax=10.0, duration=0.5, min_eclipses=3):
        times, fluxes, ferrs = self.get_lightcurve()
        if times is None:
            raise ValueError(f'no fits data found for TIC {self.tess_id:010d}.')

        # run BLS:
        periodogram, stats = run_bls(times, fluxes, ferrs, pmin, duration)

        # find the time of the first eclipse:
        # bjd0 = times[0] + period*(times[0]-times[0]//period)

        # find the log-likelihood for the first eclipse:
        # bls_logp = logfaps[max_power_idx]

        # instantiate the BLS_run instance:
        # bls_run = BLS_run(
        #     tic=self,
        #     duration=duration,
        #     min_eclipses=min_eclipses,
        #     pmin=pmin, pmax=pmax,
        #     bjd0=bjd0, period=period,
        #     max_power=max_power,
        #     bls_logp=bls_logp)
        # bls_run.save()

        # return bls_run

    @classmethod
    def add_to_triage(cls, tess_id, ephemeris, source):
        """
        Ingest an entry into triage. Thus class method does the following steps:
        * looks up TIC instance from the passed `tess_id`; if not found, it adds it via MAST;
        * TODO: checks data availability against MAST and downloads any missing data;
        * takes proposed ephemerides as a list of quadruplets (bjd0, period, model, version).

        Parameters
        ----------
        * tess_id
        * ephemerides

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
        """

        # look up tess_id and add it if not found:
        tics = cls.objects.filter(tess_id=tess_id)
        if len(tics) == 0:
            tic = cls.from_mast(tess_id)
            tic.save()
        else:
            tic = tics[0]

        # look up ephemeris source and add it if not found:
        if isinstance(source, dict):
            model = source.setdefault('model', '')
            version = source.setdefault('version', '')
            author = source.setdefault('author', '')
            reference = source.setdefault('reference', '')

            matches = EphemerisSource.objects.filter(model=model, version=version)
            if len(matches) == 0:
                source = EphemerisSource(model=model, version=version, author=author, reference=reference)
                source.save()
            elif len(matches) > 1:
                raise ValueError('source is not uniquely identified by the passed dictionary fields.')
            else:
                source = matches[0]
        elif isinstance(source, EphemerisSource):
            # we got passed the instance of EphemerisSource already, nothing to do.
            pass
        else:
            raise ValueError('source needs to be either an EphemerisSource instance or a dictionary.')

        # instantiate ephemeris if necessary:
        if isinstance(ephemeris, dict):
            bjd0 = ephemeris.setdefault('bjd0', 0.0)
            bjd0_uncert = ephemeris.setdefault('bjd0_uncert', 0.0)
            period = ephemeris.setdefault('period', 1.0)
            period_uncert = ephemeris.setdefault('period_uncert', 0.0)

            ephemeris = Ephemeris(tic=tic, bjd0=bjd0, bjd0_uncert=bjd0_uncert, period=period, period_uncert=period_uncert)
            ephemeris.save()
            ephemeris.source.add(source)
        elif isinstance(ephemeris, Ephemeris):
            ephemeris.save()
            ephemeris.source.add(source)
        else:
            raise ValueError('ephemeris needs to be either an Ephemeris instance or a dictionary')

        # add ephemeris to the TIC:
        tic.ephemerides.add(ephemeris)

        # look up origin if given and add it if not found:
        # origin = kwargs.get('origin', None)
        # if origin is not None:
        #     origins = Origin.objects.filter(name=origin)
        #     if len(origins) == 0:
        #         origin = Origin(name=origin)
        #         origin.save()
        #     else:
        #         origin = origins[0]

        # prepare the TIC for triage:
        tic.prep_for_triage(ephemerides=[ephemeris])

        return tic

    def create_static_files(self, static_dir='static/catalog', provenance=None, plot_lc=True, plot_spd=True, force_overwrite=False):
        """
        Create static files for the TIC. Static files include:
        * lightcurve plots (png, in `lc_figs` directory)
        * spectral power density plots (png, in `spd_figs` directory)

        Parameters
        ----------
        static_dir : str, optional, default 'static/catalog'
            The root directory where the static files will be stored.
        provenance : str, optional, default None
            The provenance name to select the data from. If None, all
            provenances will be used.
        plot_lc : bool, optional, default True
            If True, the method will generate lightcurve plots.
        plot_spd : bool, optional, default True
            If True, the method will generate spectral power density plots.
        force_overwrite : bool, optional, default False
            If True, the method will overwrite existing files.

        Raises
        ------
        ValueError
            If no fits data are found for the TIC.
        """

        provenances = [provenance] if provenance is not None else [p.name for p in self.provenances.all()]

        if plot_lc:
            for provenance in provenances:
                fname = f'{static_dir}/lc_figs/tic{self.tess_id:010d}.{provenance}.lc.png'
                if os.path.exists(fname) and not force_overwrite:
                    continue

                data = self.get_lightcurve(provenance=provenance)
                if len(data) <= 1:
                    continue

                plt.figure('lcfig', figsize=(16, 5))
                plt.xlabel('Truncated Barycentric Julian Date')
                plt.ylabel('Normalized PDC flux')
                plt.plot(data['times'], data['fluxes'], 'b.')
                plt.savefig(fname)
                plt.close()

                flt = data['times'] < data['times'][0] + 10  # first 10 days of data
                plt.figure('zlcfig', figsize=(16, 5))
                plt.xlabel('Truncated Barycentric Julian Date')
                plt.ylabel('Normalized PDC flux')
                plt.plot(data['times'][flt], data['fluxes'][flt], 'b.')
                plt.savefig(f'{static_dir}/lc_figs/tic{self.tess_id:010d}.{provenance}.zlc.png')
                plt.close()

        if plot_spd:
            for provenance in provenances:
                fname = f'{static_dir}/spd_figs/tic{self.tess_id:010d}.{provenance}.spd.png'
                if os.path.exists(fname) and not force_overwrite:
                    continue

                data = self.get_spd(provenance=provenance)
                if np.all(np.isnan(data[provenance]['powers'])):
                    continue

                plt.figure('spfig', figsize=(8, 5))
                plt.yscale('log')
                plt.xlabel('Period [d]')
                plt.ylabel('Lomb-Scargle power (log scale)')
                plt.plot(data[provenance]['periods'], data[provenance]['powers'], 'b-')
                plt.savefig(f'{static_dir}/spd_figs/tic{self.tess_id:010d}.{provenance}.spd.png')
                plt.close()

    def prep_for_triage(self, static_dir='static/catalog', filelist=None, ephemerides=None):
        # under review.
        return
        # lcfn = f'{static_dir}/lc_data/tic{self.tess_id:010d}.norm.lc'
        # if os.path.exists(lcfn):
        #     times, fluxes = np.loadtxt(lcfn, usecols=(0, 1), unpack=True)
        # else:
        #     times, fluxes, _ = read_from_all_fits(self.tess_id, filelist=filelist)

        # if times is None or fluxes is None:
        #     raise ValueError(f'cannot find any local data associated with TIC {self.tess_id:010d}.')

        # if ephemerides is None:
        #     ephemerides = self.ephemerides.all()

        # for eph in ephemerides:
        #     for period_str, period_mult, in zip(['half', 'period', 'double'], [0.5, 1.0, 2.0]):
        #         phases = bjd2phase(times=times, bjd0=eph.bjd0 if eph.bjd0 is not None else 0.0, period=eph.period*period_mult)

        #         plt.figure()
        #         plt.xlabel('Phase')
        #         plt.ylabel('Normalized PDC flux')
        #         plt.title(f'TIC {self.tess_id:010d}, period {eph.period*period_mult:0.8f} days')
        #         plt.plot(phases, fluxes, 'b.')
        #         plt.plot(phases[phases > +0.4]-1.0, fluxes[phases > +0.4], 'c.')
        #         plt.plot(phases[phases < -0.4]+1.0, fluxes[phases < -0.4], 'c.')
        #         plt.savefig(f'{static_dir}/triage/tic{self.tess_id:010d}.{eph.pk}.{period_str}.phase.png')
        #         plt.close()

    def __str__(self):
        return '%10d' % self.tess_id

    def __repr__(self):
        return (
            f"<TIC(tess_id={self.tess_id}, "
            f"ra={self.ra}, "
            f"dec={self.dec}, "
            f"glon={self.glon}, "
            f"glat={self.glat}, "
            f"Tmag={self.Tmag}, "
            f"teff={self.teff}, "
            f"logg={self.logg}, "
            f"abun={self.abun}, "
            f"pmra={self.pmra}, "
            f"pmdec={self.pmdec}, "
            f"gaia_id={self.gaia_id}, "
            f"gaia_dr3_id={self.gaia_dr3_id}, "
            f"kepler_id={self.kepler_id}, "
            f"asas_id={self.asas_id}, "
            f"sectors={self.sectors.all()}, "
            f"provenances={self.provenances.all()}, "
            f"datatype={self.datatype})>"
        )


class BLS_run(models.Model):
    # inputs:
    tic = models.ForeignKey('catalog.TIC', on_delete=models.PROTECT)
    duration = models.FloatField('duration', null=False, blank=False)
    min_eclipses = models.IntegerField('minimum number of eclipses', null=False, blank=False)
    pmin = models.FloatField('minimum period', null=False, blank=False)
    pmax = models.FloatField('maximum period', null=False, blank=False)

    # main deliverables:
    bjd0 = models.FloatField('time of first eclipse', null=False, blank=False)
    period = models.FloatField('period', null=False, blank=False)
    max_power = models.FloatField('maximum SPD power', null=False, blank=False)
    bls_logp = models.FloatField('log-likelihood for the first eclipse', null=False, blank=False)

    # diagnostics:
    date_run = models.DateTimeField('BLS run date', auto_now_add=True)
    depth = models.FloatField('eclipse depth', null=False, blank=False)
    depth_unc = models.FloatField('eclipse depth uncertainty', null=False, blank=False)
    depth_odd = models.FloatField('double-period odd eclipse depth', null=False, blank=False)
    depth_odd_unc = models.FloatField('double-period odd eclipse depth uncertainty', null=False, blank=False)
    depth_even = models.FloatField('double-period even eclipse depth', null=False, blank=False)
    depth_even_unc = models.FloatField('double-period even eclipse depth uncertainty', null=False, blank=False)
    depth_half = models.FloatField('half-period eclipse depth', null=False, blank=False)
    depth_half_unc = models.FloatField('half-period eclipse depth uncertainty', null=False, blank=False)

    def __str__(self):
        return f'TIC {self.tic.tess_id} t0={self.bjd0:6.6f} P0={self.period:6.6f}'

    def __repr__(self):
        return (
            f"<BLS_run(tic={self.tic}, "
            f"duration={self.duration}, "
            f"min_eclipses={self.min_eclipses}, "
            f"pmin={self.pmin}, "
            f"pmax={self.pmax}, "
            f"bjd0={self.bjd0}, "
            f"period={self.period}, "
            f"max_power={self.max_power}, "
            f"bls_logp={self.bls_logp}, "
            f"date_run={self.date_run}, "
            f"depth={self.depth}, "
            f"depth_unc={self.depth_unc}, "
            f"depth_odd={self.depth_odd}, "
            f"depth_odd_unc={self.depth_odd_unc}, "
            f"depth_even={self.depth_even}, "
            f"depth_even_unc={self.depth_even_unc}, "
            f"depth_half={self.depth_half}, "
            f"depth_half_unc={self.depth_half_unc})>"
        )


class Origin(models.Model):
    name = models.CharField(max_length=32)

    def __str__(self):
        return self.name

    def __repr__(self):
        return (
            f"<Origin(name='{self.name}')>"
        )


class EB(models.Model):
    tic = models.ForeignKey('catalog.TIC', on_delete=models.PROTECT)
    origin = models.ManyToManyField(Origin)
    sectors = models.ManyToManyField(Sector, blank=True)
    signal_id = models.IntegerField('signal id', default=1)
    in_catalog = models.BooleanField('in catalog')
    date_added = models.DateTimeField('date added', auto_now_add=True)
    date_modified = models.DateTimeField('date modified', null=True, blank=True, auto_now=True)

    # define data source choices:
    class ObjectSource(models.TextChoices):
        FFI = 'FFI', 'Full Frame Image (FFI)'
        TPF = 'TPF', 'Target Pixel File (TPF)'
        LCF = 'LCF', 'Light Curve File (LCF)'
    source = models.CharField(max_length=3, choices=ObjectSource.choices, default=ObjectSource.LCF)

    ephemeris = models.OneToOneField('catalog.Ephemeris', on_delete=models.SET_NULL, null=True, blank=True, related_name='ephemeris')
    morph_coeff = models.FloatField('morphology coefficient', null=True, blank=True)

    # eclipse properties:
    # TODO: create a new class for eclipse properties!
    prim_width_pf = models.FloatField('primary eclipse width from the polyfit model', null=True, blank=True, editable=False)
    sec_width_pf = models.FloatField('secondary eclipse width from the polyfit model', null=True, blank=True, editable=False)
    prim_depth_pf = models.FloatField('primary eclipse depth from the polyfit model', null=True, blank=True, editable=False)
    sec_depth_pf = models.FloatField('secondary eclipse depth from the polyfit model', null=True, blank=True, editable=False)
    prim_pos_pf = models.FloatField('primary eclipse position from the polyfit model', null=True, blank=True, editable=False)
    sec_pos_pf = models.FloatField('secondary eclipse position from the polyfit model', null=True, blank=True, editable=False)

    prim_width_2g = models.FloatField('primary eclipse width from the 2-gaussian model', null=True, blank=True, editable=False)
    sec_width_2g = models.FloatField('secondary eclipse width from the 2-gaussian model', null=True, blank=True, editable=False)
    prim_depth_2g = models.FloatField('primary eclipse depth from the 2-gaussian model', null=True, blank=True, editable=False)
    sec_depth_2g = models.FloatField('secondary eclipse depth from the 2-gaussian model', null=True, blank=True, editable=False)
    prim_pos_2g = models.FloatField('primary eclipse position from the 2-gaussian model', null=True, blank=True, editable=False)
    sec_pos_2g = models.FloatField('secondary eclipse position from the 2-gaussian model', null=True, blank=True, editable=False)

    # flags:
    ambiguous = models.BooleanField('ambiguous', default=False)
    insufficient = models.BooleanField('insufficient cycles', default=False)
    heartbeat = models.BooleanField('heartbeat star', default=False)
    pulscomp = models.BooleanField('pulsating component', default=False)
    multi = models.BooleanField('multiperiodic system', default=False)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['tic', 'signal_id'], name='eb_id')
        ]
        verbose_name_plural = 'EBs'

    def save(self, *args, **kwargs):
        # we are overloading the save() method to be able to test for
        # duplicates, handle many-to-many relationships, to generate
        # supporting plots and to auto-increment the id.

        # if there are no args/kwargs, revert to built-in save():
        if len(args) == 0 and len(kwargs) == 0:
            super(EB, self).save()
            return

        origin = kwargs.get('origin', None)
        tol = kwargs.get('tol', 1e-3)
        sectors = kwargs.get('sectors', [])

        force_adding = kwargs.get('force_adding', False)

        # If the TIC does not exist in the database, save it
        # irrespective of whether the period is passed or not:
        if EB.objects.filter(tic=self.tic).count() == 0:
            super(EB, self).save()
            self.origin.add(origin)
            for sector in sectors:
                self.sectors.add(sector)
            return

        # If period is not provided, and the TIC exists in the
        # database, refuse to add another one unless forced:
        if self.period is None and EB.objects.filter(tic=self.tic).count() != 0:
            if not force_adding:
                print('TIC %010d is already present in the database, not saving (pass force_adding=True to override).' % self.tic.tess_id)
                return
            else:
                self.signal_id = EB.objects.filter(tic=self.tic).count() + 1
                super(EB, self).save()
                self.origin.add(origin)
                for sector in sectors:
                    self.sectors.add(sector)

        # If period is provided and the TIC does exist in the
        # database, check its value against the tolerance:
        if self.period is not None and EB.objects.filter(tic=self.tic).count() > 0:
            counter = 1
            for entry in EB.objects.filter(tic=self.tic):
                # if the existing entry does not have a period, update the record:
                if entry.period is None:
                    break
                if abs(self.period-entry.period) < tol and not force_adding:
                    print('TIC %010d with period %5.5f-d already exists (signal %02d), not saving (pass force_adding=True to override).' % (self.tic.tess_id, self.period, counter))
                    return
                counter += 1

            self.signal_id = counter
            super(EB, self).save()
            self.origin.add(origin)
            for sector in sectors:
                self.sectors.add(sector)
            return

    def create_static_files(self, static_dir='static/catalog', provenance=None, plot_ph=True, force_overwrite=False):
        provenances = [provenance] if provenance is not None else [p.name for p in self.tic.provenances.all()]

        if plot_ph:
            bjd0 = self.bjd0 if self.bjd0 is not None else 0.0
            period = self.period if self.period is not None else 1.0

            for provenance in provenances:
                fname = f'{static_dir}/lc_figs/tic{self.tic.tess_id:010d}.{self.signal_id:02d}.{provenance}.ph.png'
                if os.path.exists(fname) and not force_overwrite:
                    continue

                data = self.tic.get_lightcurve(provenance=provenance)
                if len(data) <= 1:
                    continue

                phases = backend.bjd2phase(times=data['times'], bjd0=bjd0, period=period)

                plt.figure('phfig', figsize=(8, 5))
                plt.xlabel('Phase')
                plt.ylabel('Normalized PDC flux')
                plt.title(f'TIC {self.tic.tess_id:010d}, period {self.period:0.4f} days')
                plt.plot(phases, data['fluxes'], 'b.')
                plt.plot(phases[phases > +0.4]-1.0, data['fluxes'][phases > +0.4], 'b.')
                plt.plot(phases[phases < -0.4]+1.0, data['fluxes'][phases < -0.4], 'b.')
                plt.savefig(fname)
                plt.close()

    def __str__(self):
        return '%010d:%02d' % (self.tic.tess_id, self.signal_id)

    def __repr__(self):
        return (
            f"<EB(tic={self.tic}, "
            f"signal_id={self.signal_id}, "
            f"in_catalog={self.in_catalog}, "
            f"date_added={self.date_added}, "
            f"date_modified={self.date_modified}, "
            f"source='{self.source}', "
            f"ephemeris={self.ephemeris}, "
            f"morph_coeff={self.morph_coeff}, "
            f"prim_width_pf={self.prim_width_pf}, "
            f"sec_width_pf={self.sec_width_pf}, "
            f"prim_depth_pf={self.prim_depth_pf}, "
            f"sec_depth_pf={self.sec_depth_pf}, "
            f"prim_pos_pf={self.prim_pos_pf}, "
            f"sec_pos_pf={self.sec_pos_pf}, "
            f"prim_width_2g={self.prim_width_2g}, "
            f"sec_width_2g={self.sec_width_2g}, "
            f"prim_depth_2g={self.prim_depth_2g}, "
            f"sec_depth_2g={self.sec_depth_2g}, "f""
            f"prim_pos_2g={self.prim_pos_2g}, "
            f"sec_pos_2g={self.sec_pos_2g}, "
            f"ambiguous={self.ambiguous}, "
            f"insufficient={self.insufficient}, "
            f"heartbeat={self.heartbeat}, "
            f"pulscomp={self.pulscomp}, "
            f"multi={self.multi})>"
        )


class EphemerisSource(models.Model):
    """
    Ephemeris source model.

    The model is used to store information about the source of the ephemerides.

    Attributes
    ----------
    model : str
        The model used to generate the ephemerides.
    version : str
        The version of the model.
    author : str
        The author of the model.
    reference : str
        The reference (publication) to the model.
    """

    model = models.CharField(max_length=16)
    version = models.CharField(max_length=16, null=True, blank=True)
    author = models.CharField(max_length=64)
    reference = models.CharField(max_length=64, null=True, blank=True)

    def __str__(self):
        if self.model == 'manual':
            return f'{self.author} ({self.model})'
        return f'{self.model} {self.version}'

    def __repr__(self):
        return (
            f"<EphemerisSource(model='{self.model}', "
            f"version='{self.version}', "
            f"author='{self.author}', "
            f"reference='{self.reference}')>"
        )


class Ephemeris(models.Model):
    date_added = models.DateTimeField('date added', auto_now_add=True)
    source = models.ForeignKey(EphemerisSource, on_delete=models.CASCADE)

    tic = models.ForeignKey(TIC, related_name='ephemerides', on_delete=models.CASCADE)
    # eb = models.ForeignKey(EB, related_name='ephemerides', on_delete=models.CASCADE)

    bjd0 = models.FloatField('bjd0', null=True, blank=True)
    bjd0_uncert = models.FloatField('bjd0 uncertainty', null=True, blank=True)
    period = models.FloatField('orbital period', null=True, blank=True)
    period_uncert = models.FloatField('period uncertainty', null=True, blank=True)

    triage_timestamp = models.DateTimeField('triaged on', null=True, blank=True)
    triage_status = models.CharField(max_length=16, null=False, blank=False, default='triage')
    triage_period = models.CharField(max_length=16, null=False, blank=False, default='period')
    triage_username = models.CharField(max_length=32, null=True, blank=True)

    class Meta:
        verbose_name_plural = 'Ephemerides'

    def __str__(self):
        return '%s + E x %s' % (self.bjd0, self.period)

    def __repr__(self):
        return (
            f"<Ephemeris(date_added={self.date_added}, "
            f"source={self.source}, "
            f"tic={self.tic}, "
            f"bjd0={self.bjd0}, "
            f"bjd0_uncert={self.bjd0_uncert}, "
            f"period={self.period}, "
            f"period_uncert={self.period_uncert}, "
            f"triage_timestamp={self.triage_timestamp}, "
            f"triage_status='{self.triage_status}', "
            f"triage_period='{self.triage_period}', "
            f"triage_username='{self.triage_username}')>"
        )

    def generate_plots(self, static_dir='static/catalog', provenance=None, force_overwrite=False):
        """
        Generate plots for the ephemeris.

        Parameters
        ----------
        static_dir : str, optional, default 'static/catalog'
            The root directory where the static files will be stored.
        provenance : str, optional, default None
            The provenance name to select the data from. If None, all
            provenances will be used.
        force_overwrite : bool, optional, default False
            If True, the method will overwrite existing files.
        """

        bjd0 = self.bjd0 or 0.0
        period = self.period or 1.0

        data = backend.load_data(self.tic.tess_id, provenance=provenance)
        if len(data) <= 1:
            raise ValueError(f'no lightcurve data found for TIC {self.tic.tess_id:010d}.')

        fname = f'{static_dir}/triage/tic{self.tic.tess_id:010d}.{provenance or "all"}.eph{self.pk:07d}.ph.png'
        if os.path.exists(fname) and not force_overwrite:
            return

        backend.plot_plc(
            filename=fname,
            data=data,
            bjd0=bjd0,
            period=period,
            title=f'TIC {self.tic.tess_id:010d}, period {self.period:0.4f} days',
            xlabel='Phase',
            ylabel='Normalized PDC flux'
        )


class Comment(models.Model):
    author = models.CharField('author', null=False, blank=False, max_length=32)
    text = models.TextField('text', null=True, blank=True)
    timestamp = models.DateTimeField('commented on')

    ephem = models.ForeignKey(Ephemeris, related_name='comments', null=True, on_delete=models.CASCADE)
    eb = models.ForeignKey(EB, related_name='comments', null=True, on_delete=models.CASCADE)
    tic = models.ForeignKey(TIC, related_name='comments', null=True, on_delete=models.CASCADE)

    def __str__(self):
        return "%s's comment on %s" % (self.author, self.timestamp)

    def __repr__(self):
        return (
            f"<Comment(author='{self.author}', "
            f"text='{self.text}', "
            f"timestamp={self.timestamp}, "
            f"ephem={self.ephem}, "
            f"eb={self.eb}, "
            f"tic={self.tic})>"
        )
