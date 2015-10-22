from __future__ import absolute_import

import pickle

from scipy.io import loadmat
from numpy import zeros, where

from django.forms import Form, ModelForm, CharField, TextInput, ChoiceField, ValidationError
from django.core.exceptions import ValidationError

#from django.db.models import Model, ForeignKey
from django.core.urlresolvers import reverse
from django.template.defaultfilters import slugify

from crispy_forms.helper import FormHelper
from crispy_forms.bootstrap import StrictButton, FieldWithButtons, AppendedText, TabHolder, Tab
from crispy_forms.layout import Layout, ButtonHolder, Submit, Div, Field, Fieldset, HTML

#from selectable.forms import AutoCompleteSelectField
#from .lookups import GeneSymbolLookup

from sapling.settings import BASE_DIR

from .models import Experiment, GFP


class ExperimentForm(ModelForm):
  class Meta:
    fields = ('Annotation', 'Name',)
    model = Experiment
	
  def __init__(self, *args, **kwargs):
    super(ExperimentForm, self).__init__(*args, **kwargs)
    self.fields['Annotation'].label = "Gene list: insert as a list of gene symbol names (one gene per line)"
    self.fields['Name'].label = "Name of the new experiment"
    self.helper = FormHelper()
    self.helper.layout = Layout('Annotation', 'Name',	
      ButtonHolder(Submit('create', 'Create', css_class='btn-primary')),)

  def clean_Annotation(self):
    mapper = BASE_DIR + "/data/Mappers/Mapper_Symbol_Entrez.pkl"
    handle = open(mapper, 'r')
    symbolsentrez = pickle.load(handle)
    handle.close()

    # ProteinID
    path = BASE_DIR + "/data/entrezGeneID_Human_1_24_2015.mat"
    proteinid = loadmat(path, squeeze_me=False, mat_dtype=True)
    proteinid = proteinid['proteinid']

    data = self.cleaned_data['Annotation']

    genes = data.split()

    # Annotation
    annotation = zeros(proteinid.shape)
    missingmap = []
    missingnet = []
    for protein in genes:
      if len(protein) == 1:
        print protein
      if protein in symbolsentrez.keys():
        if int(symbolsentrez[protein]) in proteinid:
          if len(protein) == 1:
            print protein
          index = where(int(symbolsentrez[protein]) == proteinid)[0]
          annotation[index] = 1
        else:
          missingnet.append(proteinid)
      else:
        missingmap.append(protein)	

    if sum(annotation) < 20:
      raise ValidationError("Your gene list has less than 20 genes in the network")
    elif sum(annotation) > 300:
      raise ValidationError("Yur genes list has more than 300 genes in the network")
    return data

  def clean_Name(self):
    data = self.cleaned_data['Name']
    print len(Experiment.objects.filter(Name=self.cleaned_data['Name']))
    if len(Experiment.objects.filter(Name=self.cleaned_data['Name'])) > 0:
      raise ValidationError("This experiment name is already used")
    return data


class ExperimentRunForm(ModelForm):
  class Meta:
    #fields = ('Name',)
    exclude = ('Name', 'Annotation',)
    model = Experiment

  def __init__(self, *args, **kwargs):
    super(ExperimentRunForm, self).__init__(*args, **kwargs)
    self.helper = FormHelper()
    self.helper.layout = Layout( 
      ButtonHolder(Submit('run', 'Run', css_class='btn_primary')),)


class GFPForm(ModelForm):
  networkid = CharField(max_length=512, required=False)
  PPI_CHOICES = (
    ('none', '-----------'),
    ('BioGRID', 'BioGRID'),
    ('HIPPIE', 'HIPPIE'),
    ('IntAct', 'IntAct'),
    ('I2D', 'I2D'),
    ('GM', 'GeneMANIA'))
  ppi = ChoiceField(choices=PPI_CHOICES, required=False)
  SM_CHOICES = (
    ('none', '-----------'),
    ('KEGG', 'KEGG'),
    ('Reactome', 'Reactome'),
    ('Phenocarta', 'Phenocarta'),
    ('InterPro', 'InterPro'),
    ('Pfam', 'Pfam'))
  sm = ChoiceField(choices=SM_CHOICES, required=False)

  class Meta:
    fields = ('Algorithm', 'networkid', 'ppi', 'sm')
    model = GFP

  def __init__(self, *args, **kwargs):
    super(GFPForm, self).__init__(*args, **kwargs)
    self.fields['networkid'].label = "Network ID"
    self.fields['ppi'].label = "Protein-protein interaction database"
    self.fields['sm'].label = "Semantic similarity database"
    self.helper = FormHelper()
    self.helper.form_class = 'blueForms'
    self.helper.layout = Layout( \
      Fieldset("Function inference algorithm", 'Algorithm'), \
      Fieldset("Navigate network", \
        TabHolder(Tab("Protein-protein interaction", 'ppi'), \
          	Tab("Co-expression", HTML(" <p><a href={% url 'gfps:experiments:annotation' %} class='btn'>Find a network ID</a></p>"), 'networkid'), \
		Tab("Semantic similarity", 'sm'))), \
		ButtonHolder(Submit('add', 'Add', css_class='btn-primary')))

