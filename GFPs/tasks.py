from __future__ import absolute_import

import logging
from datetime import datetime

from django.template.defaultfilters import slugify
from django.core.mail import EmailMessage
from django.utils import timezone
from django.utils.timezone import activate

from djcelery import celery

from algorithms.algofact import Algorithm
from report.report import Report
from sapling.settings import BASE_DIR, TIME_ZONE
from .models import ProcessManagement, Experiment, GFP, Performance

log = logging.getLogger(__name__)

activate(TIME_ZONE)

REPORT_DIR = BASE_DIR + "/report/"

@celery.task(bind=True, name='gfppredict')
def predict(self, gfpID, user):
  """
    annotation: list of gene symbols
  """
  log.debug("Start time: %s" % timezone.now().strftime('%d/%m/%Y %H:%M:%S'))

  gfp = GFP.objects.filter(id=gfpID)
  gfp = gfp[0]

  log.debug("RUNNING: initialize method")
  self.update_state(state='RUNNING', meta='Initialize method')
  if not ProcessManagement.objects.filter(gfp_id=gfpID, experiment_id=gfp.experiment_id).exists():
    ProcessManagement.objects.create(gfp_id=gfpID, experiment_id=gfp.experiment_id, \
      State= "RUNNING: initialize method", Name=gfp.Network+"_"+gfp.Algorithm)
  else:
    ProcessManagement.objects.filter(gfp_id=gfpID, experiment_id=gfp.experiment_id).update(State= \
      "RUNNING: initialize method")

  algorithm = Algorithm()
  algorithm = algorithm.getAlgorithm(gfp.Algorithm)

  log.debug("RUNNING: load files")
  self.update_state(state='RUNNING', meta='Load files')
  ProcessManagement.objects.filter(gfp_id=gfpID, experiment_id=gfp.experiment_id).update(State= \
    "RUNNING: load files")

  algorithm.loadData(gfpID)

  log.debug("RUNNING: run method")
  self.update_state(state='RUNNING', meta='Run method')
  ProcessManagement.objects.filter(gfp_id=gfpID, experiment_id=gfp.experiment_id).update(State= \
    "RUNNING: run method")
  algorithm.run()

  log.debug("RUNNING: save predictions")
  self.update_state(state='RUNNING', meta='Save predictions')
  ProcessManagement.objects.filter(gfp_id=gfpID, experiment_id=gfp.experiment_id).update( \
    State="RUNNING: save predictions")

  algorithm.save(gfpID)

  log.debug("RUNNING: enrichment")
  algorithm.enrichment(gfp.experiment_id, gfpID)

  ProcessManagement.objects.filter(gfp_id=gfpID, experiment_id=gfp.experiment_id).update( \
    State="FINISHED")
  log.debug("FINISHED")

  gfps = GFP.objects.filter(experiment_id=gfp.experiment_id)

  processes = ProcessManagement.objects.filter(experiment_id=gfp.experiment_id)

  flag = True
  if len(processes) != len(gfps):
    flag = False
  else:
    for process in processes:
      if process.State != "FINISHED":
        flag = False
        break

  log.debug("Experiment flag: %s", str(flag))

  ## Create report and aggregate
  if flag:
    log.debug("REPORT EMAIL")

    gfps = GFP.objects.filter(experiment_id=gfp.experiment_id)

    # Aggregation of GFPs
    algorithm.aggregate(gfps, gfp.experiment_id)

    experiment = Experiment.objects.filter(id=gfp.experiment_id)
    log.debug("RUNNING: start report")
    report = Report(gfp.experiment_id)
    filename = REPORT_DIR + "sapling_gfp_report_" + slugify(experiment[0].Name) + "_" \
		+ str(gfp.experiment_id) + ".pdf"
    log.debug("RUNNING: %s" % filename)
    report.generateReport(filename)
    del report

    subject = "SAPLING report: " + experiment[0].Name
    body = "Dear " + user.first_name + " "+ user.last_name + ",\n\n"
    body += "Please find attached a report of your experiment " + experiment[0].Name + ".\n\n"
    body += "Regards,\n SAPLING"
	
    email = EmailMessage(subject=subject, body=body, from_email="sapling.gfp@gmail.com", to=[user.email])
    reportname = "sapling_gfp_report_" + slugify(experiment[0].Name) + "_" + str(gfp.experiment_id) + ".pdf"
    handle = open(filename, 'rb')
    email.attach(reportname, handle.read(), "application/pdf")

    enrichname = "Enrichment_GeneList_GO_" + str(gfp.experiment_id) + ".csv"
    filename = REPORT_DIR + "Enrichment_GeneList_GO_" + str(gfp.experiment_id) + ".csv"
    handle = open(filename, 'rb')
    email.attach(enrichname, handle.read(), "text/csv")

    email.send()
    handle.close()
    log.debug("DONE")


@celery.task(bind=True, name='aggregate')
def aggregate(self, experimentID, user):
  """
    Aggregate the GFP scores and send report
  """
  gfps = GFP.objects.filter(experiment_id=experimentID)
  log.debug("RUNNING aggregate: initialize method")

  algorithm = Algorithm()
  algorithm = algorithm.getAlgorithm("NV")
  algorithm.loadDataAggregate()
  algorithm.aggregate(gfps, experimentID)

  experiment = Experiment.objects.filter(id=experimentID)
  log.debug("RUNNING aggregate: start report")
  report = Report(experimentID)
  filename = REPORT_DIR + "sapling_gfp_report_" + slugify(experiment[0].Name) + "_" + str(experimentID) + ".pdf"
  report.generateReport(filename)
  del report

  subject = "SAPLING report: " + experiment[0].Name
  body = "Dear " + user.first_name + " "+ user.last_name + ",\n\n"
  body += "Please find attached a report of your experiment " + experiment[0].Name + ".\n\n"
  body += "Regards,\n SAPLING"

  email = EmailMessage(subject=subject, body=body, from_email="sapling.gfp@gmail.com", to=[user.email])
  reportname = "sapling_gfp_report_" + slugify(experiment[0].Name) + "_" + str(experimentID) + ".pdf"
  handle = open(filename, 'rb')
  email.attach(reportname, handle.read(), "application/pdf")

  enrichname = "Enrichment_GeneList_GO_" + str(experimentID) + ".csv"
  filename = REPORT_DIR + "Enrichment_GeneList_GO_" + str(experimentID) + ".csv"
  handle = open(filename, 'rb')
  email.attach(enrichname, handle.read(), "text/csv")

  email.send()
  handle.close()

  log.debug("DONE Aggregate")
