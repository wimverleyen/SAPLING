from sys import exc_info

from django.contrib.auth.models import User
from django.contrib import messages
from django.db.models import Model, CharField, TextField, SlugField, IntegerField, FloatField, DateTimeField, ImageField, ForeignKey
from django.db.utils import IntegrityError
from django.core.urlresolvers import reverse
from django.template.defaultfilters import slugify


class Experiment(Model):
  """
    Experiment: list of GFP methods
  """
  user = ForeignKey(User, related_name='experiments')
  Annotation = TextField()
  Name = CharField(max_length=255)
  slug = SlugField(max_length=255, blank=True)

  def save(self, *args, **kwargs):
    self.slug = slugify(self.Name)
    super(Experiment, self).save(*args, **kwargs)
  def get_absolute_url(self):
    return reverse('gfps:experiments:detail', kwargs={'slug': self.slug})
  def get_absolute_url_run(self):
    return reverse('gfps:experiments:run', kwargs={'slug': self.slug})
  def running_gfp(self):
    process = ProcessManagement.objects.filter(experiment_id=self, State__icontains='RUNNING')
    return len(process)
  def active_processes(self):
    active = 0
    try:
      dj = djcelery.celery
      i = dj.control.inspect()
      active = len(i.active().items()[0][1])
      log.debug('Experiment Run detail view get context data active processes: %d', len(i.active().items()[0][1]))
    except Exception:
      log.error('Experiment Run detail view: not able proceesses')
    return active


class GeneSymbol(Model):
  Symbol = CharField(max_length=255)
  ProteinID = IntegerField()


class CoAnnotation(Model):
  Geoid = CharField(max_length=255)
  Category = CharField(max_length=255)
  Description = CharField(max_length=1023)
  Samples = IntegerField()
  Type = CharField(max_length=255)


class GFP(Model):
  """
    GFP method definition
  """
  # (Database, User)
  ALGORITHM_CHOICES = (
    ('NV', 'Neighbor Voting'),
    ('LR', 'Logistic Regression'),
    ('RWR', 'Random Walk with Restart'),
    ('SGD', 'Stochastic Gradient Descent'),
    ('PA', 'Passive Aggressive'),
    ('RP', 'RankProp'))

  TYPE_CHOICES = (('PPI', 'Protein-protein interaction'),
    ('Co - RNASeq', 'Co-expression: RNA-Seq'),
    ('Co - Microarray', 'Co-expression: Microarray'),
    ('SM', 'Semantic similarity'))	
	
  experiment = ForeignKey(Experiment, related_name='gfps')
  Network = CharField(max_length=255)
  Name = CharField(max_length=255)
  Type = CharField(max_length=255, choices=TYPE_CHOICES)	
  Algorithm = CharField(max_length=255, choices=ALGORITHM_CHOICES)
  slug = SlugField(max_length=255, blank=True)

  class Meta:
    unique_together = ('experiment', 'slug')

  def __unicode__(self):
    return "GFP: " + self.Algorithm + "_" + self.Network

  def save(self, request, form, *args, **kwargs):
    if form.data['ppi'] != 'none':
      #print "PPI network"
      self.Network = form.data['ppi']
      self.Type = "PPI"
    elif form.data['networkid'] != '':
      #print "Co network"
      self.Network = form.data['networkid']
      self.Type = "Co"
    elif form.data['sm'] != 'none':
      self.Network = form.data['sm']
      self.Type = "SM"

    self.Name = "GFP: " + self.Algorithm + "_" + self.Network
    self.slug = slugify(self.Algorithm+self.Network)

    try:
      super(GFP, self).save(*args, **kwargs)
    except IntegrityError as e:
      print "Exception Integrity: ", exc_info()[0]
      print "Exception Integrity: ", e
      messages.warning(request, "Please do not use the same algorithm and data resource again")
    except:
      print "Exception: ", exc_info()[0]

  def store(self, experiment, network, type, algorithm, *args, **kwargs):
    self.experiment = experiment
    self.Network = network
    self.Type = type
    self.Algorithm = algorithm
    self.Name = "GFP: " + self.Algorithm + "_" + self.Network
    self.slug = slugify(self.Algorithm+self.Network)

    try:
      super(GFP, self).save(*args, **kwargs)
    except IntegrityError as e:
      print "Exception Integrity: ", exc_info()[0]
      print "Exception Integrity: ", e
      messages.warning(request, "Please do not use the same algorithm and data resource again")
    except:
      print "Exception: ", exc_info()[0]

  def get_absolute_url_prediction(self):
    return reverse('gfps:experiments:prediction', kwargs={'slug1': self.experiment.slug, \
					'slug2':self.slug})


class Performance(Model):
  gfp = ForeignKey(GFP, related_name='performance')
  CV = IntegerField()
  AUROC = FloatField()

  def save(self, gfpID, auroc, *args, **kwargs):
    self.gfp_id = gfpID
    self.AUROC = auroc
    self.CV = 3
    super(Performance, self).save(*args, **kwargs)

  def get_experiment(self):
    return self.gfp.experiment.Name


class Prediction(Model):
  gfp = ForeignKey(GFP, related_name='predictor')
  CV = IntegerField()
  ProteinID = IntegerField()
  Score = FloatField()
  Rank= FloatField()
  Annotation = IntegerField()

  def save(self, gfpID, ProteinID, score, rank, annotation, CV=-1, *args, **kwargs):
    self.gfp_id = gfpID
    self.CV=CV
    self.ProteinID=ProteinID
    self.Score=score
    self.Rank=rank
    self.Annotation=annotation
    super(Prediction, self).save(*args, **kwargs)


class NovelCandidates(Model):
  gfp = ForeignKey(GFP, related_name='novelcandidates')
  CV = IntegerField()
  Symbol = CharField(max_length=255)
  ProteinID = IntegerField()
  Score = FloatField()
  Rank = FloatField()
  Annotation = IntegerField()


class Enrichment(Model):
  experiment = ForeignKey(Experiment, related_name='enrichexperiment')
  Profile = CharField(max_length=255)
  ProfileID = IntegerField()
  Pvalue = FloatField()
  BHPvalue = FloatField(null=True)


class EnrichmentGFP(Model):
  gfp = ForeignKey(GFP, related_name='enrichgfp')
  Profile = CharField(max_length=255)
  ProfileID = IntegerField()
  Test = CharField(max_length=255)
  Pvalue = FloatField()
  BHPvalue = FloatField(null=True)


class Homologene(Model):
  ProteinIDHuman = IntegerField()
  SymbolHuman = CharField(max_length=255)
  ProteinIDMouse = IntegerField()
  SymbolMouse = CharField(max_length=255)


class GO(Model):
  GOID = IntegerField()
  Name = CharField(max_length=255)
  Namespace = CharField(max_length=255)


class KEGG(Model):
  KEGGID = IntegerField()
  Name = CharField(max_length=255)


class Phenocarta(Model):
  CartaID = IntegerField()
  Name = CharField(max_length=255)
  ExternalID = CharField(max_length=255, null=True)


class AggregatePerformance(Model):
  experiment = ForeignKey(Experiment, related_name='aggregateperformance')
  CV = IntegerField()
  AUROC = FloatField()

  def save(self, experimentID, auroc, *args, **kwargs):
    self.experiment_id = experimentID
    self.AUROC = auroc
    self.CV = 3
    super(AggregatePerformance, self).save(*args, **kwargs)

  def get_experiment(self):
    return self.gfp.experiment.name


class AggregatePrediction(Model):
  experiment = ForeignKey(Experiment, related_name='aggregatepredictor')
  CV = IntegerField()
  ProteinID = IntegerField()
  Score = FloatField()
  Rank = FloatField()
  Annotation = IntegerField()

  def save(self, experimentID, ProteinID, score, rank, annotation, CV=-1, *args, **kwargs):
    self.experiment_id = experimentID
    self.CV=CV
    self.ProteinID=ProteinID
    self.Score=score
    self.Rank=rank
    self.Annotation=annotation
    super(AggregatePrediction, self).save(*args, **kwargs)


class CMA(Model):
  experiment = ForeignKey(Experiment, related_name='cma')
  Configuration = CharField(max_length=511)
  Performance = FloatField()
  Reproducibility = FloatField()
  Similarity = FloatField(default=0.0)


class ProcessManagement(Model):
  gfp = ForeignKey(GFP, related_name='process')
  experiment = ForeignKey(Experiment, related_name='experiment')
  Name = CharField(max_length=255)
  State = CharField(max_length=255)
  Time = DateTimeField(blank=True, null=True)
