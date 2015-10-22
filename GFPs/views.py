from __future__ import absolute_import

from datetime import datetime
import logging

from django.db.models import Count
from django.contrib import messages
from django.shortcuts import redirect
from django.http import HttpResponse, Http404
from django.views.generic import View, ListView, CreateView, DetailView, RedirectView
from django.contrib.auth.models import User

import djcelery

from braces.views import LoginRequiredMixin, SetHeadlineMixin, PrefetchRelatedMixin

from .models import CoAnnotation, Experiment, GFP, Performance, Prediction, ProcessManagement
from .forms import ExperimentForm, GFPForm, ExperimentRunForm
from .tasks import predict, aggregate
from report.report import Report

log = logging.getLogger(__name__)

class RestrictToUserMixin(LoginRequiredMixin):
  def get_queryset(self):
    queryset = super(RestrictToUserMixin, self).get_queryset()
    queryset = queryset.filter(user=self.request.user)
    return queryset


class AnnotationListView(ListView):
  model = CoAnnotation
  headline = 'Annotations'

  def get_queryset(self):
    queryset = super(AnnotationListView, self).get_queryset()
    return queryset


class ExperimentView(View):
  def get(self, request, *args, **kwargs):
    return HttpResponse("An experiment view")


class ExperimentListView(RestrictToUserMixin, LoginRequiredMixin, ListView):
  model = Experiment
  headline = 'Experiments'
  prefetch_related = ('gfps', 'process')

  def get_queryset(self):
    log.debug('Experiment list view get queryset')
    queryset = super(ExperimentListView, self).get_queryset()
    queryset = queryset.annotate(gfp_count=Count('gfps'))
    return queryset
  

class ExperimentCreateView(RestrictToUserMixin, LoginRequiredMixin, SetHeadlineMixin, CreateView):
  form_class = ExperimentForm
  headline = 'Create'
  model = Experiment

  def form_valid(self, form):
    self.object = form.save(commit=False)
    self.object.user = self.request.user
    if not Experiment.objects.filter(Name=form.data['Name']).exists():
      self.object.save()
      log.debug('Experiment create view form_valid: %s created', form.data['Name'])
      return super(ExperimentCreateView, self).form_valid(form)
    else:
      messages.warning(self.request, "Your experiment already exists")
      log.debug('Experiment create view form_valid: %s already exists', form.data['Name'])
      return redirect('gfps:experiments:list')


class ExperimentDetailView(RestrictToUserMixin, LoginRequiredMixin, PrefetchRelatedMixin, DetailView):
  form_class = GFPForm
  http_method_names = ['get', 'post']
  model = Experiment
  prefetch_related = ('gfps',)

  def get_context_data(self, **kwargs):
    log.debug('Experiment detail view get')
    context = super(ExperimentDetailView, self).get_context_data(**kwargs)
    context.update({'form': self.form_class(self.request.POST or None)})
    return context
	
  def post(self, request, *args, **kwargs):
    log.debug('Experiment detail view post')
    form = self.form_class(request.POST)
    if form.is_valid():
      log.debug('Experiment detail view post: valid form')
      log.debug('Experiment detail view post: valid form PPI: %s', form.data['ppi'])
      log.debug('Experiment detail view post: valid form SM: %s', form.data['sm'])
      log.debug('Experiment detail view post: valid form networkid: %s', form.data['networkid'])
      obj = self.get_object()
      if form.data['ppi'] != 'none':
        log.debug('Experiment detail view post: PPI network')
        gfp = form.save(commit=False)
        gfp.experiment = obj
        gfp.save(request, form)
      elif form.data['sm'] != 'none':
        log.debug('Experiment detail view post: SM network')
        gfp = form.save(commit=False)
        gfp.experiment = obj
        gfp.save(request, form)
      elif form.data['networkid'] != '':
        log.debug('Experiment detail view post: Co network')
        ids = CoAnnotation.objects.filter(Geoid = form.data['networkid'])
        if len(ids) > 0:
          id = ids[0]
          log.debug('Experiment detail view post: good network ID')
          log.debug('Experiment detail view post: GEOID: %s' % id.Geoid)
          gfp = form.save(commit=False)
          gfp.experiment = obj
          gfp.save(request, form)
        else:
          log.debug('Experiment detail view post: wrong network ID')
          messages.warning(request, "Please use another co-expression network identifier")
    else:
      log.debug('Experiment detail view post: form not valid')
      return self.get(request, *args, **kwargs)
    return redirect(obj)


class ExperimentRunDetailView(RestrictToUserMixin, LoginRequiredMixin, PrefetchRelatedMixin, DetailView):
  form_class = ExperimentRunForm
  http_method_names = ['get', 'post']
  model = Experiment
  prefetch_related = ('gfps',)
  template_name = "GFPdevs/experiment_run.html"

  def get_context_data(self, **kwargs):
    log.debug('Experiment Run detail view get context data')
    context = super(ExperimentRunDetailView, self).get_context_data(**kwargs)
    try:
      dj = djcelery.celery
      i = dj.control.inspect()
      log.debug('Experiment Run detail view get context data active processes: %d', len(i.active().items()[0][1]))
    except Exception:
      log.error('Experiment Run detail view: not able proceesses')

    context.update({'form': self.form_class(self.request.POST or None)})
    return context
	
  def post(self, request, *args, **kwargs):
    log.debug('Experiment Run detail view post')

    #Click on run bottom
    obj = self.get_object()
	
    log.debug('Experiment Run detail view post experiment name: %s' % obj.Name)

    block = False
    processes = ProcessManagement.objects.filter(experiment_id=obj.id)
    log.debug('Experiment Run detail view post %d in process management' % len(processes))

    for process in processes:
      if "RUNNING" in process.State:
        block = True
        break

    if block:
      print "RUN BLOCK"
      log.debug('Experiment Run detail view post BLOCK RUNNING')
      messages.warning(request, "Please wait until your results are emailed to you")
      return redirect('gfps:experiments:list')
    else:
      flag = True
      log.debug('Experiment Run detail view post %d GFPs' % len(obj.gfps.all()))

      for gfp in obj.gfps.all():
        log.debug('Experiment Run detail view post %d id GFPs' % gfp.id)

        #if not Performance.objects.filter(gfp=gfp).exists():
        if not ProcessManagement.objects.filter(gfp_id=gfp.id).exists():
          result = predict.delay(gfp.id, self.request.user)
          log.debug('Experiment Run detail view post task id: %s' % result.task_id)
          flag = False
      if flag and len(obj.gfps.all()) > 0:
        log.debug('Experiment Run detail view post aggregate')
        aggregate.delay(obj.id, self.request.user)
      return redirect('gfps:experiments:list')


class ExperimentRemoveGFPView(LoginRequiredMixin, RedirectView):
  model = GFP

  def get_redirect_url(self, *args, **kwargs):
    log.debug('Experiment remove GFP redirect URL')
    return self.experiment.get_absolute_url()
		
  def get_object(self, pk, experiment_pk):
    try:
      log.debug('Experiment remove GFP get object')
      gfp = self.model.objects.get(pk=pk, experiment_id=experiment_pk, \
					experiment__user=self.request.user)
    except GFP.DoesNotExist:
      log.error('Experiment remove GFP not found')
      raise Http404
    else:
      log.error('Experiment remove GFP return')
      return gfp

  def get(self, request, *args, **kwargs):
    self.object = self.get_object(kwargs.get('pk'), kwargs.get('experiment_pk'))
    log.debug('Experiment remove GFP get')
    self.experiment = self.object.experiment
    messages.success(request, u'{0.Name} was removed from {1.Name}'.format(self.object, self.experiment))
    self.object.delete()
    return super(ExperimentRemoveGFPView, self).get(request, *args, **kwargs)


class PerformanceListView(RestrictToUserMixin, LoginRequiredMixin, PrefetchRelatedMixin, ListView):
  model = Performance
  headline = 'Performance'
  prefetch_related = ('experiment', 'gfp')
  template_name = "GFPs/experiment_performance.html"

  def get_queryset(self):
    experiment = Experiment.objects.filter(slug=self.kwargs['slug'], user=self.request.user)
    gfps = GFP.objects.filter(experiment=experiment)
    return Performance.objects.filter(gfp=gfps)


class PredictionListView(RestrictToUserMixin, LoginRequiredMixin, PrefetchRelatedMixin, ListView):
  model = Performance
  headline = 'Prediction'
  prefetch_related = ('performance', )
  template_name = "GFPs/experiment_prediction.html"

  def get_queryset(self):
    experiment = Experiment.objects.filter(slug=self.kwargs['slug1'], user=self.request.user)
    gfp = GFP.objects.filter(experiment=experiment, slug=self.kwargs['slug2']) 
    return Prediction.objects.filter(gfp=gfp, Annotation=0).order_by('Score')[0:50]
