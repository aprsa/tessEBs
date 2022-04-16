from django.urls import path, include
from django.views.generic import TemplateView, DetailView

from . import views

urlpatterns = [
	path('', views.MainView.as_view(), name='main'),
    path('<int:tess_id>.<int:signal_id>', views.DetailsView.as_view(), name='details'),
    path('<int:tess_id>', views.DetailsView.as_view(), name='details'),
    path('about', TemplateView.as_view(template_name='static/about.html'), name='about'),
    path('help', TemplateView.as_view(template_name='static/help.html'), name='help'),
    path('mast', views.MASTExport.as_view(), name='mast_export'),
    path('support', views.SupportPDF.as_view(), name='letter_of_support'),
    path('search', views.SearchView.as_view(), name='search'),
    path('search_results', views.SearchResultsView.as_view(), name='search_results'),
    path('triage', views.react_triage),
    path('triage/<int:tess_id>', views.react_triage),
    path('triage_done', TemplateView.as_view(template_name='static/triage_done.html'), name='triage_done'),
    #path('triage/user', views.react_triage),
    path('api/triage/next', views.api_triage_next),
    path('api/triage/ephem/<int:tess_id>/<str:username>', views.api_triage_ephem),
    path('api/triage/ephem/<int:tess_id>', views.api_triage_ephem),
    path('api/triage/save', views.api_triage_save),
    path('api/triage/user/<str:username>', views.api_triage_user),
    path('api/data/lc/<int:tess_id>', views.api_data_lc),
    path('api/data/periodogram/<int:tess_id>', views.api_data_periodogram),
    path('accounts/', include('django.contrib.auth.urls')),
    path('accounts/profile/', views.profile, name='profile'),
]
