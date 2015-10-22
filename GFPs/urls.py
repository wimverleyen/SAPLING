from __future__ import absolute_import

from django.conf.urls import patterns, url, include

from .views import ExperimentView, ExperimentListView, ExperimentCreateView, ExperimentDetailView,\
			ExperimentRunDetailView, ExperimentRemoveGFPView, PerformanceListView, \
			PredictionListView, AnnotationListView

experiment_patterns = patterns(
	'',
	url(
		regex=r'^$', 
		view = ExperimentListView.as_view(), 
		name='list'),
	url(
		regex=r'^annotation/$', 
		view = AnnotationListView.as_view(), 
		name='annotation'),
	url(
		regex=r'^create/$', 
		view = ExperimentCreateView.as_view(), 
		name='create'),
	url(
		regex=r'^d/(?P<slug>[-\w]+)/$', 
		view = ExperimentDetailView.as_view(), 
		name='detail'),
	#url(
	#	regex=r'^n/(?P<slug>[-\w]+)/$', 
	#	view = NetworkDetailView.as_view(), 
	#	name='net'),
	url(
		regex=r'^r/(?P<slug>[-\w]+)/$', 
		view = ExperimentRunDetailView.as_view(), 
		name='run'),
	url(
		#regex=r'^pe/(?P<experiment_pk>[-\w]+)/$', 
		regex=r'^pe/(?P<slug>[-\w]+)/$', 
		view = PerformanceListView.as_view(), 
		name='performance'),
	url(
		regex=r'^pr/(?P<slug1>[-\w]+)/(?P<slug2>[-\w]+)/$', 
		view = PredictionListView.as_view(), 
		name='prediction'),
	url(
		regex=r'^remove/(?P<experiment_pk>\d+)/(?P<pk>\d+)$',
		view=ExperimentRemoveGFPView.as_view(),
		name='remove_gfp'),
)

urlpatterns = patterns(
	'',
	url(r'experiments/', include(experiment_patterns, namespace='experiments')),
)
