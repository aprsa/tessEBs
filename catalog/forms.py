from django import forms
from django.utils.safestring import mark_safe

class ChangeForm(forms.Form):
    is_eb = forms.BooleanField(label='Is EB', required=False)
    bjd0 = forms.FloatField(label='BJD0', required=False)
    period = forms.FloatField(label='Period', required=False)

class SearchForm(forms.Form):
    order_options = [
        ('',            ''),
        ('period',      'Period'),
        ('morph_coeff', 'Morphology coefficient'),
        ('tic__ra',     'Right ascension'),
        ('tic__dec',    'Declination'),
        ('tic__glon',   'Galactic longitude'),
        ('tic__glat',   'Galactic latitude'),
        ('tic__teff',   'Temperature'),
        ('tic__logg',   'Surface gravity'),
        ('tic__abun',   'Chemical abundance'),
        ('tic__Tmag',   'TESS magnitude')
    ]

    order_by = forms.ChoiceField(label='Order by:', required=True, initial='period', choices=order_options)
    then_order_by = forms.ChoiceField(label='then order by:', required=False, initial='', choices=order_options)

    period_min = forms.FloatField(label=mark_safe('P<sub>0</sub> between'), required=False, min_value=0, max_value=10000)
    period_max = forms.FloatField(label='', required=False, min_value=0, max_value=10000)