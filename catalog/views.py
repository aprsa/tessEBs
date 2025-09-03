import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re
import csv
import json

import numpy as np
from django.conf import settings
from django.shortcuts import render, redirect
from django.http import JsonResponse, FileResponse
from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required
from django.views.decorators.csrf import csrf_exempt
from django.core.paginator import Paginator
from django.db.models import Q, Count
from django.contrib.auth import logout
from django.urls import reverse_lazy

from csv_export.views import CSVExportView

from django.views.generic import View, TemplateView, ListView
from django.views.generic.detail import BaseDetailView
from django.contrib.auth.views import LogoutView, PasswordChangeView, PasswordResetView
from django.contrib.auth.forms import PasswordChangeForm, PasswordResetForm

from . import backend
from . import pipeline as pl
from .ephemeris_ui import create_ephemeris_ui
from .models import EB, TIC, Ephemeris, EphemerisSource, Comment
from .forms import ChangeForm

import json
from scipy.signal import find_peaks
from random import randint
from datetime import datetime


class MainView(ListView):
    template_name = 'catalog/mainlist.html'

    def __init__(self):
        EBs = EB.objects.filter(in_catalog=True).order_by('tic__tess_id')
        self.ebno = EBs.count()
        self.paginator = Paginator(EBs, 100)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['link'] = '?order_by=tic'
        context['ebno'] = self.ebno
        context['page_of_EBs'] = self.paginator.get_page(self.page_num)
        context['fields'] = ['tic__tess_id', 'sectors', 'morph_coeff', 'source', 'flags']
        return context

    def get_queryset(self):
        self.page_num = self.request.GET.get('page', 1)


class DetailsView(TemplateView):
    # Yes, it would be tempting to have this subclassed from DetailView, but that isn't easy.
    # The problem is that we don't use pk as an identifier and it gets messier than simply
    # using the TemplateView.

    model = EB
    template_name = 'catalog/details.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        tess_id = context.get('tess_id', None)
        signal_id = context.get('signal_id', 1)
        context['signal_id'] = signal_id

        try:
            tic = TIC.objects.get(tess_id=tess_id)
        except TIC.DoesNotExist:
            tic = None

        context['tic'] = tic
        context['list_of_ebs'] = EB.objects.filter(tic__tess_id=tess_id)

        if tic is not None:
            # check if the syndicated fits file exists:
            fits_file = os.path.join(settings.BASE_DIR, 'static', 'catalog', 'lc_data', f'tic{tess_id:010d}.fits')
            context['data_exist'] = os.path.exists(fits_file)
        else:
            context['data_exist'] = False

        return context


class SupportPDF(BaseDetailView):
    def get_context_data(self, **kwargs):
        context = {}
        return context

    def get(self, request, *args, **kwargs):
        context = self.get_context_data()
        response = FileResponse(open('static/support.pdf', 'rb'), content_type='application/pdf')
        response['Content-Disposition'] = "inline; filename='static/support.pdf'"
        return response


class MASTExport(CSVExportView):
    model = EB
    verbose_names = False
    specify_separator = False
    filename = 'mast_data_export.csv'
    fields = (
        'tic__tess_id',
        'signal_id',
        'date_added',
        'date_modified',
        'source',
        'tic__ra',
        'tic__dec',
        'tic__pmra',
        'tic__pmdec',
        'tic__Tmag',
        'tic__teff',
        'ephemeris__bjd0',
        'ephemeris__bjd0_uncert',
        'ephemeris__period',
        'ephemeris__period_uncert',
        'morph_coeff',
        'prim_width_pf',
        'prim_depth_pf',
        'prim_pos_pf',
        'sec_width_pf',
        'sec_depth_pf',
        'sec_pos_pf',
        'prim_width_2g',
        'prim_depth_2g',
        'prim_pos_2g',
        'sec_width_2g',
        'sec_depth_2g',
        'sec_pos_2g',
        'sectors',
    )

    def get_queryset(self):
        return self.model.objects.filter(in_catalog=True)

    def get_field_value(self, obj, field_name):
        # overloading this to strip timezone info.
        if field_name == 'date_modified' or field_name == 'date_added':
            field = obj._meta.get_field(field_name)
            val = field.value_from_object(obj)
            if val:
                val = val.replace(tzinfo=None)
            return val
        return super().get_field_value(obj, field_name)

    def get_csv_writer_fmtparams(self):
        fmtparams = super().get_csv_writer_fmtparams()
        fmtparams['dialect'] = 'unix'
        fmtparams['quoting'] = csv.QUOTE_MINIMAL
        return fmtparams


class SearchView(TemplateView):
    template_name = 'catalog/search.html'


class SearchResultsView(ListView):
    model = EB
    template_name = 'catalog/search_results.html'

    def filter(self):
        self.fields = [
            'tic__tess_id',
            'ephemeris__bjd0',
            'ephemeris__period',
            'morph_coeff',
            'tic__ra',
            'tic__dec',
            'tic__glon',
            'tic__glat',
            'tic__teff',
            'tic__logg',
            'tic__abun',
            'tic__Tmag',
            'sectors__sector_id',
            'nsectors',
            'signal_id'
        ]

        filter_list = {}
        for field in self.fields:
            ll = self.request.GET.get(field+'_min', None)
            ul = self.request.GET.get(field+'_max', None)
            if ll is not None and ll != '':
                filter_list[f'{field}__gte'] = float(ll)
            if ul is not None and ul != '':
                filter_list[f'{field}__lte'] = float(ul)
        return filter_list

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['format'] = self.request.GET.get('display_format', 'html')
        context['ebno'] = self.ebno
        if context['format'] == 'html':
            context['page_of_EBs'] = self.page_of_EBs
        else:
            context['page_of_EBs'] = self.object_list
        context['link'] = re.sub('&page=[0-9]+', '', self.search_criteria)
        context['fields'] = self.request.GET.getlist('c')
        # if quicksearch is used, no fields will be selected; fall back:
        if len(context['fields']) == 0:
            context['fields'] = ['tic__tess_id', 'sectors', 'ephemeris__bjd0', 'ephemeris__bjd0_uncert', 'ephemeris__period', 'ephemeris__period_uncert', 'morph_coeff', 'source', 'flags']

        return context

    def get_queryset(self):
        # store search criteria for pagination:
        self.search_criteria = self.request.get_full_path()
        # self.search_criteria = self.request.GET.dict()
        # self.search_criteria = {c: search_criteria[c][0] for c in search_criteria.keys()}

        # check if quicksearch (by TIC) called the view:
        tic_query = self.request.GET.get('tic', None)
        if tic_query:
            self.object_list = EB.objects.filter(Q(tic__tess_id__icontains=tic_query))

        else:
            order_by = self.request.GET.get('order_by')
            then_order_by = self.request.GET.get('then_order_by')
            ordering = [order_by, then_order_by] if then_order_by != 'none' else [order_by]

            filter_list = self.filter()

            incat_only = self.request.GET.get('incat_only')
            if incat_only == 'on':
                filter_list['in_catalog'] = True

            self.object_list = EB.objects.annotate(nsectors=Count('sectors')).filter(**filter_list).order_by(*ordering)

        self.ebno = self.object_list.count()
        self.paginator = Paginator(self.object_list, 100)
        page_num = self.request.GET.get('page', 1)
        self.page_of_EBs = self.paginator.get_page(page_num)

        return self.object_list


class EphemView(TemplateView):
    template_name = 'catalog/ephem.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        tess_id = context.get('tess_id', None)
        signal_id = context.get('signal_id', 1)

        script, div = create_ephemeris_ui(tess_id)

        context['signal_id'] = signal_id
        tic_list = TIC.objects.filter(tess_id=tess_id)
        if tic_list.count() == 0:
            context['tic'] = None
        else:
            context['tic'] = tic_list[0]
        context['list_of_ebs'] = EB.objects.filter(tic__tess_id=tess_id)

        context['bokeh_script'] = script
        context['bokeh_div'] = div

        return context


class TriageView(TemplateView):
    template_name = 'catalog/triage.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        tess_id = context.get('tess_id', None)

        tic = TIC.objects.get(tess_id=tess_id)
        context['tic'] = tic
        context['list_of_ebs'] = EB.objects.filter(tic__tess_id=tess_id)

        # check if the syndicated fits file exists:
        fits_file = os.path.join(settings.BASE_DIR, 'static', 'catalog', 'lc_data', f'tic{tess_id:010d}.fits')
        context['data_exist'] = os.path.exists(fits_file)

        # check if static plots exist:
        context['plots'] = {}
        for provenance in tic.provenances.all():
            context['plots'][provenance.name] = {}
            lc_file = os.path.join(settings.BASE_DIR, 'static', 'catalog', 'lc_figs', f'tic{tess_id:010d}.{provenance.name}.lc.png')
            zlc_file = os.path.join(settings.BASE_DIR, 'static', 'catalog', 'lc_figs', f'tic{tess_id:010d}.{provenance.name}.zlc.png')
            spd_file = os.path.join(settings.BASE_DIR, 'static', 'catalog', 'spd_figs', f'tic{tess_id:010d}.{provenance.name}.spd.png')
            ph_file = os.path.join(settings.BASE_DIR, 'static', 'catalog', 'lc_figs', f'tic{tess_id:010d}.{provenance.name}.ph.png')

            context['plots'][provenance.name]['lc'] = f'catalog/lc_figs/{os.path.basename(lc_file)}' if os.path.exists(lc_file) else None
            context['plots'][provenance.name]['zlc'] = f'catalog/lc_figs/{os.path.basename(zlc_file)}' if os.path.exists(zlc_file) else None
            context['plots'][provenance.name]['spd'] = f'catalog/spd_figs/{os.path.basename(spd_file)}' if os.path.exists(spd_file) else None
            context['plots'][provenance.name]['ph'] = f'catalog/lc_figs/{os.path.basename(ph_file)}' if os.path.exists(ph_file) else None

        return context


class ChangePasswordView(PasswordChangeView):
    form_class = PasswordChangeForm
    success_url = reverse_lazy('profile')
    template_name = 'registration/change_password.html'


class ResetPasswordView(PasswordResetView):
    form_class = PasswordResetForm
    success_url = reverse_lazy('login')
    template_name = 'registration/reset_password.html'


class GoodbyeView(LogoutView):
    template_name = 'registration/goodbye.html'

    def post(self, request):
        logout(request)
        return render(request, 'registration/goodbye.html')


class TestApiView(TemplateView):
    template_name = 'catalog/test_api.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        return context


@method_decorator(login_required, name='dispatch')
class ApiTestView(View):
    """
    Tests the API.

    Returns:
    --------
    JsonResponse
    """

    def post(self, request):
        return JsonResponse({
            'status': 'success',
        })


@method_decorator(login_required, name='dispatch')
class ApiTestErrorView(View):
    """
    Tests the API error handling.
    """
    
    def post(self, request):
        raise ValueError("This is a test error for the API.")


@method_decorator(login_required, name='dispatch')
class ApiTicDownloadMetaView(View):
    """
    Downloads metadata for a given TIC ID from MAST and SIMBAD.
    """

    def post(self, request):
        data = request.POST

        tess_id = data.get('tess_id', None)
        if tess_id is None:
            return JsonResponse({
                'status': 'failure',
                'error': 'tess_id not provided'
            }, status=400)

        meta = backend.download_meta(int(tess_id))

        return JsonResponse({
            'status': 'success',
            **meta
        })


@method_decorator(login_required, name='dispatch')
class ApiSyndicateDataView(View):
    """
    Syndicate data for a given TIC ID.

    POST parameters:
    ----------------
    tess_id: int
        TIC ID to syndicate data for.
    force_overwrite: bool
        If True, overwrite existing data.

    Returns:
    --------
    JsonResponse
    """

    def post(self, request):
        data = request.POST

        tess_id = data.get('tess_id', None)
        if tess_id is None:
            return JsonResponse({
                'status': 'failure',
                'error': 'tess_id not provided'
            }, status=400)

        try:
            tic = TIC.objects.get(tess_id=tess_id)
        except TIC.DoesNotExist:
            return JsonResponse({
                'status': 'failure',
                'error': f'tess_id {tess_id} not found'
            }, status=400)

        force_overwrite = data.get('force_overwrite', False)
        tic.syndicate_data(force_overwrite=force_overwrite)
        tic.attach_spds_to_data()

        return redirect(request.META.get('HTTP_REFERER'))
        return JsonResponse({
            'status': 'success',
        })


@method_decorator(login_required, name='dispatch')
class ApiLombScargleView(View):
    def post(self, request):
        try:
            data = json.loads(request.body)

            time = np.array(data.get('time'))
            flux = np.array(data.get('flux'))
            pmin = data.get('pmin', 0.1)
            pmax = data.get('pmax', 10.0)
            pstep = data.get('pstep', 0.001)
            npeaks = data.get('npeaks', 0)

            data = {
                'times': time,
                'fluxes': flux
            }
            # np.savetxt('data.txt', data)
            pgram = pl.run_lombscargle(data, pmin, pmax, pstep, npeaks)

            return JsonResponse({
                'period': pgram['periods'].tolist(),
                'power': pgram['powers'].tolist()
            })
        except Exception as e:
            return JsonResponse({'error': str(e)}, status=400)


@method_decorator(login_required, name='dispatch')
class ApiTicAddView(View):
    def post(self, request):
        try:
            data = request.POST

            tess_id = data.get('tess_id', None)
            if tess_id is None:
                return JsonResponse({
                    'status': 'failure',
                    'error': 'tess_id not provided'
                }, status=400)

            TIC.from_mast(int(tess_id), syndicate_data=False, create_static=False)
            return redirect(request.META.get('HTTP_REFERER'))
            # comment out above and uncomment below to debug
            # return JsonResponse({
            #     'status': 'success',
            #     'data': data,
            #     'tess_id': tess_id,
            #     # 'data': data,
            # })
        except Exception as e:
            return JsonResponse({
                'status': 'failure',
                'error': str(e),
            }, status=400)


@method_decorator(csrf_exempt, name='dispatch')
class ApiEphemAddView(View):
    def post(self, request):
        try:
            data = json.loads(request.body.decode('utf-8'))

            tess_id = data.get('tess_id')
            t0 = data.get('t0')
            period = data.get('period')
            model = data.get('model')
            version = data.get('version')
            reference = data.get('reference')
            attribution = data.get('attribution')

            if model == 'manual' and attribution == 'self':
                attribution = f'{request.user.first_name} {request.user.last_name}'

            eph_source, created = EphemerisSource.objects.get_or_create(author=attribution, model=model, version=version, reference=reference)
            if created:
                eph_source.save()

            tic = TIC.objects.get(tess_id=tess_id)
            ephem = Ephemeris(tic=tic, source=eph_source, bjd0=t0, period=period)
            ephem.save()

            if data.get('generate_plots', True):
                ephem.generate_plots()

            return JsonResponse({
                'status': 'success',
            })
        except Exception as e:
            return JsonResponse({
                'status': 'failure',
                'error': str(e)
            }, status=400)







def profile(request):
    context = {}
    return render(request, 'registration/profile.html', context)


def get_changes(request):
    if request.method == 'POST':
        form = ChangeForm(request.POST)
        if form.is_valid():
            # do something
            return HttpResponseRedirect('/thanks/')
    else:
        form = ChangeForm()


def _get_random_ephem_triage():
    # TODO: filter by triage status instead of all
    ephems_triage = Ephemeris.objects.filter(triage_status='triage', tic__datatype='LCF')
    if ephems_triage.count() <= 1:
        return None

    return ephems_triage[randint(0, ephems_triage.count()-1)]


########## REACT VIEWS #############

def react_triage(request, tess_id=None):
    if tess_id is None:
        triage_entry = _get_random_ephem_triage()
        if triage_entry is None:
            return redirect('/triage_done')
        # tess_id = triage_entry.eb.tic.tess_id
        tess_id = triage_entry.tic.tess_id
        return redirect('/triage/{}'.format(tess_id))
    if request.user.is_authenticated:
        return render(request, "react/index.html", {'username': request.user.username})
    else:
        return redirect('/admin/login/?next=/triage/{}'.format(tess_id))


########## API VIEWS ###############

def api_triage_next(request):
    return JsonResponse({'tic': _get_random_ephem_triage().eb.tic.tess_id})


def _api_phase_up(t, t0=0, P=1):
    return (t-t0) % P / P


def _api_ephem(tess_id, username=None):
    def _ephem_dict(ephem, i):
        source = ephem.source.all()[0]
        return {
            'pk': ephem.pk,
            'i': i,
            'source': source.model if len(source.model) else source.author,
            'period': ephem.period,
            'bjd0': ephem.bjd0,
            'triage_status': ephem.triage_status,
            'triage_period': ephem.triage_period,
            'comments': [{'author': comment.author, 'text': comment.text} for comment in ephem.comments.all()]
        }

    tic = TIC.objects.get(tess_id=tess_id)
    ephems = tic.ephemerides.all()

    # eb = EB.objects.get(tic__tess_id=tess_id, signal_id=1)
    # ephems = eb.ephemerides.all()

    all_ephems = Ephemeris.objects.all()
    ephem_triage_total = all_ephems.count()
    ephem_triage_done = all_ephems.exclude(triage_status='triage').count()
    if username:
        ephem_triage_user = all_ephems.filter(triage_username=username).exclude(triage_status='triage').count()
    else:
        ephem_triage_user = 0

    return {
        'tic': tess_id,
        'ephems': [_ephem_dict(ephem, i+1) for i, ephem in enumerate(ephems)],
        'ephem_triage_total': ephem_triage_total,
        'ephem_triage_done': ephem_triage_done,
        'ephem_triage_user': ephem_triage_user
    }


def api_triage_ephem(request, tess_id, username=None):
    return JsonResponse(_api_ephem(tess_id, username))


def api_triage_user(request, username):
    ephems = Ephemeris.objects.filter(triage_username=username).order_by('-triage_timestamp')
    triaged = []
    tics = []
    for e in ephems:
        # if e.eb.tic.tess_id not in tics:
        #     tics.append(e.eb.tic.tess_id)
        #     triaged.append({'tic': e.eb.tic.tess_id, 'timestamp': e.triage_timestamp.strftime('%Y-%m-%d %H:%M')})
        if e.tic.tess_id not in tics:
            tics.append(e.tic.tess_id)
            triaged.append({'tic': e.tic.tess_id, 'timestamp': e.triage_timestamp.strftime('%Y-%m-%d %H:%M')})

    comments = []
    for c in Comment.objects.filter(author=username).order_by('-timestamp'):
        if c.ephem is not None:
            # comments.append({'tic': c.ephem.eb.tic.tess_id, 'comment': c.text, 'timestamp': c.timestamp.strftime('%Y-%m-%d %H:%M')})
            comments.append({'tic': c.ephem.tic.tess_id, 'comment': c.text, 'timestamp': c.timestamp.strftime('%Y-%m-%d %H:%M')})

    return JsonResponse({'triaged': triaged, 'comments': comments})


@csrf_exempt
def api_triage_save(request):
    if request.method != 'POST':
        return JsonResponse({'success': False, 'msg': 'only accepts POST'})

    r = json.loads(request.body.decode('utf-8'))
    username = r.get('username')
    tess_id = r.get('tic')
    changes = r.get('changes', {})
    # TODO: remove this (or manually add usernames/permissions) once done with testing
    # if username not in ['kconroy', 'aprsa', 'andrej', 'TEST-USER']:
    #    return JsonResponse({'success': False, 'msg': 'not authorized to save changes (please report this if you think you should have permissions)', 'username': username})
    if tess_id is None:
        return JsonResponse({'success': False, 'msg': 'tess_id not provided', 'username': username})

    tess_id = int(float(tess_id))
    for ephem_pk, ephem_changes in changes.items():
        if not len(ephem_changes.keys()):
            # probably shouldn't be in the changes dictionary,
            # but let's make sure we don't waste time in case there is a blank entry
            continue

        if ephem_pk == "NEW":
            es, created = EphemerisSource.objects.get_or_create(author=username, version='(triage)', model='')
            # eb = EB.objects.get(tic__tess_id=tess_id, signal_id=1)
            tic = TIC.objects.get(tess_id=tess_id)
            e = Ephemeris(tic=tic, date_added=datetime.now())
            e.save()
            e.source.set([es])
            ephem_pk = e.pk

            period = float(ephem_changes.get('period', 1.0))
            bjd0 = float(ephem_changes.get('bjd0', 0.0))

            # need to create phase plots for this new entry
            times, fluxes = np.loadtxt('./static/catalog/lc_data/tic%010d.norm.lc' % (tess_id), usecols=(0, 1), unpack=True)
            for period_str, period_mult in zip(['half', 'period', 'double'], [0.5, 1, 2]):
                phases = _api_phase_up(times, t0=bjd0, P=period*period_mult)

                plt.clf()
                plt.xlabel('Phase')
                plt.ylabel('Normalized flux')
                plt.title('TIC %010d, period %0.8f days' % (tess_id, period*period_mult))
                plt.plot(phases, fluxes, 'b.')
                plt.plot(phases[phases > 0.9]-1.0, fluxes[phases > 0.9], 'c.')
                plt.plot(phases[phases < 0.1]+1.0, fluxes[phases < 0.1], 'c.')
                plt.savefig('./static/catalog/triage/tic%010d.%d.%s.phase.png' % (tess_id, ephem_pk, period_str))

        else:
            e = Ephemeris.objects.get(pk=int(float(ephem_pk)))

        for k, v in ephem_changes.items():
            if k == 'new_comment':
                c = Comment(ephem=e, author=username, text=v, timestamp=datetime.now())
                c.save()
            else:
                setattr(e, k, v)

        e.triage_timestamp = datetime.now()
        e.triage_username = username
        e.save()

    # return the new dictionary so changes can be confirmed and new manual entries will show static plots
    j = _api_ephem(tess_id, username)
    # also return success so react knows everything went well and to not show the error message
    j['success'] = True
    j['username'] = username
    j['msg'] = 'save succesfull! (except for comments)'
    return JsonResponse(j)


def api_data_lc(request, tess_id):

    times, fluxes = np.loadtxt('catalog/static/catalog/lc_data/tic%010d.norm.lc' % (tess_id), usecols=(0, 1), unpack=True)

    return JsonResponse({'times': times.tolist(), 'fluxes': fluxes.tolist(), 'fluxes_min': fluxes.min(), 'fluxes_max': fluxes.max()})


def api_data_periodogram(request, tess_id):
    # TODO: optimize this by downsampling the original files so we can just load and serve
    freq, lsamp, lsfap = np.loadtxt('catalog/static/catalog/spd_data/tic%010d.ls.spd' % (tess_id), unpack=True)

    # distance may need to be changed if the input is downsampled, its in units of the freq step
    peaks_inds, props = find_peaks(lsamp, height=0.001, distance=10)
    peaks_amps = lsamp[peaks_inds]
    peaks_inds_sorted = peaks_inds[peaks_amps.argsort()][::-1]
    peak_periods = 1./freq[peaks_inds_sorted]
    # peak_periods = np.linspace(0.1, 10, 101).tolist()

    # downsample if file is over 80k lines
    if len(freq) > 80000:
        freq = freq[::int(len(freq)/40000)]
        lsamp = lsamp[::int(len(lsamp)/40000)]
    periods = 1./freq

    return JsonResponse({'periods': periods.tolist(), 'powers': lsamp.tolist(), 'peak_periods': peak_periods.tolist()})
