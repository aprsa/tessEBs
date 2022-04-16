from django.db import models
import numpy as np
from astroquery.mast import Catalogs as cat
# import pipeline as p


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

    class Meta:
        verbose_name_plural = 'TICs'

    @classmethod
    def from_mast(cls, tess_id):
        tic = cls(tess_id=tess_id)
        mast_entry = cat.query_criteria(catalog='TIC', ID=tess_id)
        tic.ra = mast_entry['ra'][0] if ~np.isnan(mast_entry['ra'][0]) else None
        tic.dec = mast_entry['dec'][0] if ~np.isnan(mast_entry['dec'][0]) else None
        tic.glon = mast_entry['gallong'][0] if ~np.isnan(mast_entry['gallong'][0]) else None
        tic.glat = mast_entry['gallat'][0]  if ~np.isnan(mast_entry['gallat'][0]) else None
        tic.Tmag = mast_entry['Tmag'][0] if ~np.isnan(mast_entry['Tmag'][0]) else None
        tic.teff = mast_entry['Teff'][0] if ~np.isnan(mast_entry['Teff'][0]) else None
        tic.logg = mast_entry['logg'][0] if ~np.isnan(mast_entry['logg'][0]) else None
        tic.abun = mast_entry['MH'][0] if ~np.isnan(mast_entry['MH'][0]) else None
        tic.pmra = mast_entry['pmRA'][0] if ~np.isnan(mast_entry['pmRA'][0]) else None
        tic.pmdec = mast_entry['pmDEC'][0] if ~np.isnan(mast_entry['pmDEC'][0]) else None
        tic.gaia_id = int(mast_entry['GAIA'][0]) if ~np.isnan(int(mast_entry['GAIA'][0])) else None

        return tic

    def __str__(self):
        return '%10d' % self.tess_id


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


class Origin(models.Model):
    name = models.CharField(max_length=32)

    def __str__(self):
        return self.name


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

    bjd0 = models.FloatField('bjd0', null=True, blank=True)
    bjd0_uncert = models.FloatField('bjd0 uncertainty', null=True, blank=True)
    period = models.FloatField('orbital period', null=True, blank=True)
    period_uncert = models.FloatField('period uncertainty', null=True, blank=True)
    ephem_2g_logprob = models.FloatField('2-gaussian logp', null=True, blank=True)
    ephem_pf_logprob = models.FloatField('polyfit logp', null=True, blank=True)

    morph_coeff = models.FloatField('morphology coefficient', null=True, blank=True)

    # eclipse properties:
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

#    lc = models.FilePathField('light curve', null=False, default='')

    def save(self, *args, **kwargs):
        # we are overloading the save() method to be able to test for duplicates, handle many-to-many relationships,
        # to generate supporting plots and to auto-increment the id.

        # filelist = kwargs.get('filelist', None)
        # generate_plots = kwargs.get('generate_plots', False)

        # if generate_plots:
        #     p.generate_plots(self.tic.tess_id, prefix='catalog/static/catalog/images', tlc=True, plc=True, spd=True, signal_id=self.signal_id, bjd0=self.bjd0, period=self.period, filelist=filelist)
        #     # FIXME: this shouldn't be returning
        #     return
        
        # if there are no args/kwargs, revert to built-in save():
        if len(args) == 0 and len(kwargs) == 0:
            super(EB, self).save()
            return

        origin = kwargs.get('origin', None)
        tol = kwargs.get('tol', 1e-3)
        sectors = kwargs.get('sectors', [])

        force_adding = kwargs.get('force_adding', False)

        # If the TIC does not exist in the database, save it irrespective of whether the period is passed or not:
        if EB.objects.filter(tic=self.tic).count() == 0:
            super(EB, self).save()
            self.origin.add(origin)
            for sector in sectors:
                self.sectors.add(sector)
            return

        # If period is not provided, and the TIC exists in the database, refuse to add another one unless forced:
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

        # If period is provided and the TIC does exist in the database, check its value against the tolerance:
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

    def __str__(self):
        return '%010d:%02d' % (self.tic.tess_id, self.signal_id)


class EphemerisSource(models.Model):
    model = models.CharField(max_length=16)
    version = models.CharField(max_length=16, null=True, blank=True)
    author = models.CharField(max_length=64)
    reference = models.CharField(max_length=64, null=True, blank=True)

    def __str__(self):
        return '%s %s' % (self.model, self.version)


class Ephemeris(models.Model):
    date_added = models.DateTimeField('date added', auto_now_add=True)
    source = models.ManyToManyField(EphemerisSource)

    eb = models.ForeignKey(EB, related_name='ephemerides', on_delete=models.CASCADE)

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


class Comment(models.Model):
    author = models.CharField('author', null=False, blank=False, max_length=32)
    text = models.TextField('text', null=True, blank=True)
    timestamp = models.DateTimeField('commented on')

    ephem = models.ForeignKey(Ephemeris, related_name='comments', null=True, on_delete=models.CASCADE)
    eb = models.ForeignKey(EB, related_name='comments', null=True, on_delete=models.CASCADE)
    tic = models.ForeignKey(TIC, related_name='comments', null=True, on_delete=models.CASCADE)

    def __str__(self):
        return "%s's comment on %s" % (self.author, self.timestamp)
