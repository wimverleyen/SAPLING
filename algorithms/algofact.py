from unittest import TestCase, makeSuite, main

import os
import pickle
import logging

from sklearn.linear_model import LogisticRegression, SGDClassifier, PassiveAggressiveClassifier
from math import fabs
from numpy import zeros, where, floor, asarray, dot, copy, tile, transpose, ones_like, subtract, argsort, append, delete, arange, intersect1d, linspace, spacing
from numpy.random import permutation, seed
from numpy.linalg import inv
from scipy import interp
from scipy.io import loadmat, savemat
from scipy.sparse import issparse, lil_matrix
from scipy.stats import rankdata, spearmanr
from scipy.stats.kde import gaussian_kde
from sklearn.metrics import roc_curve, auc
from matplotlib.pyplot import figure, plot, axis

from matplotlib import rc
rc('font',**{'family':'LinLibertine'})

from collections import Counter

from django.db import transaction
from sapling.settings import BASE_DIR

from enrich.enrichment import EnrichmentAnalysis

from GFPs.models import Performance, Prediction, AggregatePerformance, AggregatePrediction, GeneSymbol, NovelCandidates, GFP, Experiment, Enrichment, EnrichmentGFP, CoAnnotation

from bulk_update.helper import bulk_update

import networkx as nx
import seaborn as sns

log = logging.getLogger(__name__)

sns.set(style="white", palette="muted")

REPORT_DIR = BASE_DIR + "/report/"
 

class Algorithm:
  def __init__(self):
    seed(66)
    self._proteinid = None
    self._network = None
    self._transition = None
    self._annotation = None
    self._auroc = None
    self._scores = None

    self._foldscores = []
    self._foldlabels = []
    self._foldproteins = []

    self._fpr = None
    self._tpr = None

  def _createTransitionMatrix(self, path):
    (x, y) = self._network.shape

    D = lil_matrix(self._network.shape)

    nodedegree = self._network.sum(axis=1)
    nodedegree = asarray(nodedegree)
    nodedegree = nodedegree.ravel()
		
    D.setdiag(nodedegree)

    Q = dot(inv(D.todense()), self._network)
    print Q.shape

    savemat(path, {'Q':Q})

  def _TPR_FPR_OLD(self, annotations, scores):
    self._fpr, self._tpr, thresholds = roc_curve(annotations, scores)
    auroc = auc(self._fpr, self._tpr)

    if (self._auroc > 0.5 and auroc < 0.5) or (self._auroc < 0.5 and auroc > 0.5):
      annotations = [-1 if x==1 else x for x in annotations]
      annotations = [1 if x==0 else x for x in annotations]
      annotations = [0 if x==-1 else x for x in annotations]
		
    self._fpr, self._tpr, thresholds = roc_curve(annotations, scores)

  def _TPR_FPR(self, annotations, scores):
    scoresidx = argsort(scores)

    sortedannotations = zeros(len(annotations))
    sortedscores = zeros(len(scores))

    i = 0
    for index in scoresidx:
      sortedannotations[i] = annotations[index]
      sortedscores[i] = scores[index]
      i += 1

    tpr = [0.0]
    fpr = [0.0]

    annotations = asarray(annotations)
    annotations = annotations.ravel()
    (Npos,) = where(annotations == 1)[0].shape
    (Nneg,) = where(annotations == 0)[0].shape

    i = 0
    tp = 0
    tn = 0
    for label in sortedannotations:
      if label == 1:
        tp += 1
      else:
        tn += 1

      tpr.append(tp/float(Npos))
      fpr.append(tn/float(Nneg))

    self._fpr = fpr
    self._tpr = tpr

  def _standRank(self, scores):
    ranks = rankdata(scores)
    (x, y) = self._proteinid.shape
    ranks /= float(x)
    ranks = 1-ranks
    return ranks

  def AUROC(self, annotations, scores):
    ranks = rankdata(scores)

    Npos = len([label for label in annotations if label > 0])
    Nneg = len([label for label in annotations if label <= 0])

    ranksum = 0.0
    index = 0
    for rank in ranks:
      if annotations[index] == 1:
        ranksum += rank
      index += 1

    value = ranksum - ((Npos * (Npos + 1))/float(2))
    if Npos > 0:
      value = value/float(Npos * Nneg)
    else:
      value = 0.5
    auc = 1 - value
    return auc

  def loadData(self, gfpID, ProteinID="entrezGeneID_Human_1_24_2015.mat", mapper="Mapper_Symbol_Entrez.pkl"):

    gfp = GFP.objects.filter(id=gfpID)
    #print len(gfp)

    if len(gfp) != 1:
      log.debug("ERROR!! duplicated gfp ID")

    gfp = gfp[0]
    log.debug("GFP type= %s", gfp.Type)
    log.debug("GFP network= %s", gfp.Network)

    # Network -- REFACTORING: uniform network name
    path = ""
    if gfp.Type == 'PPI':
      path = BASE_DIR + "/data/"+ gfp.Type + "/" + gfp.Network + "/ppiNetwork_Human_1_24_2015.mat"
    elif gfp.Type == 'Co':
      path = BASE_DIR + "/data/"+ gfp.Type + "/" + gfp.Network + "/coNetwork_Human_1_24_2015.mat" 
    elif gfp.Type == 'SM':
      path = BASE_DIR + "/data/"+ gfp.Type + "/" + gfp.Network + "/smNetwork_Human_1_24_2015.mat"
    else:
      log.debug("Invalid network - type")

    network = loadmat(path, squeeze_me=False, mat_dtype=True)
    network = network['network']
    if issparse(network):
      self._network = network.todense()
    else:
      self._network = network

    path = BASE_DIR + "/data/"+ gfp.Type + "/" + gfp.Network + "/transition_Human_1_24_2015.mat"

    # transition matrix
    if os.path.exists(path):
      transition = loadmat(path, squeeze_me=False, mat_dtype=True)
      transition = transition['Q']
      if issparse(transition):
        self._transition = transition.todense()
      else:
        self._transition = transition
    else:
      self._createTransitionMatrix(path)
      transition = loadmat(path, squeeze_me=False, mat_dtype=True)
      transition = transition['Q']
      if issparse(transition):
        self._transition = transition.todense()
      else:
        self._transition = transition

    # ProteinID
    path = BASE_DIR + "/data/" + ProteinID
    proteinid = loadmat(path, squeeze_me=False, mat_dtype=True)
    self._proteinid = proteinid['proteinid']

    # Annotation
    path = BASE_DIR + "/data/Mappers/" + mapper
    handle = open(path, 'r')
    symbolsentrez = pickle.load(handle)
    handle.close()

    experiment = Experiment.objects.filter(id=gfp.experiment_id)
    experiment = experiment[0]
    annotation = experiment.Annotation.split()

    self._annotation = zeros(self._proteinid.shape)
    missingmap = []
    missingnet = []
    for protein in annotation:
      if len(protein) == 1:
        print protein
      if protein in symbolsentrez.keys():
        if int(symbolsentrez[protein]) in self._proteinid:
          if len(protein) == 1:
            print protein
          index = where(int(symbolsentrez[protein]) == self._proteinid)[0]
          self._annotation[index] = 1
          NovelCandidates.objects.get_or_create(gfp_id=gfpID, Symbol=protein, \
                            ProteinID=int(symbolsentrez[protein]), Score=0.0, \
                            Rank=0.0, Annotation=1, CV=-1)
        else:
          missingnet.append(proteinid)
      else:
        missingmap.append(protein)

  def loadDataAggregate(self, ProteinID="entrezGeneID_Human_1_24_2015.mat", mapper="Mapper_Symbol_Entrez.pkl"):
    path = BASE_DIR + "/data/" + ProteinID
    proteinid = loadmat(path, squeeze_me=False, mat_dtype=True)
    self._proteinid = proteinid['proteinid']

  def plotROCCurve(self, ID, aggregate=False, rep=[]):
    fig = figure()
    if aggregate:

      fig.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.9, hspace=0.6, wspace=0.4)
      add_plot = fig.add_subplot(1, 2, 1)

      identity = arange(0, 1.1, 0.1)

      add_plot.plot(identity, identity, color='0.5')
      add_plot.plot(self._fpr, self._tpr, color='forestgreen')

      label = "False positive rate"
      add_plot.set_xlabel(label, fontsize=12)
      label = "True positive rate"
      add_plot.set_ylabel(label, fontsize=12)
      title = "ROC curve (area = %.4f)"%self._auroc
      add_plot.set_title(title, fontsize=14)
      add_plot.set_aspect(1./add_plot.get_data_ratio())
      add_plot.grid(True)

      add_plot = fig.add_subplot(1, 2, 2)
      add_plot.set_ylim([-1.0, 1.4])
      width = 0.2
      #N = 5
      #ind = arange(N)
      #ind = [0]

      add_plot.bar(0.0, rep[0], width, color='forestgreen', label='PPI')
      add_plot.bar(0.2, rep[1], width, color='dodgerblue', label='SM')
      add_plot.bar(0.4, rep[2], width, color='gold', label='sc RNA Seq')
      add_plot.bar(0.6, rep[3], width, color='olive', label='RNA Seq')
      add_plot.bar(0.8, rep[4], width, color='indianred', label='Microarray')

      label = "Data resource"
      add_plot.set_xlabel(label, fontsize=12)
      add_plot.set_xticks([])
      label = "Total node degree\nrepresentation bias (r)"
      add_plot.set_ylabel(label, fontsize=12, multialignment='center')
      add_plot.set_aspect(1./add_plot.get_data_ratio())
      add_plot.legend(ncol=2, fontsize=6, loc="upper right")
      add_plot.grid(True)

      filename = BASE_DIR + "/figures/ROC_Curve_Aggregate_" + str(ID) + ".png"
      fig.savefig(filename, format="png", dpi=300)

    else:

      fig.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.9, hspace=0.6, wspace=0.4)
	
      add_plot = fig.add_subplot(1, 2, 1)

      identity = arange(0, 1.1, 0.1)
      add_plot.plot(identity, identity, color='0.5')
      add_plot.plot(self._fpr, self._tpr, color='forestgreen')

      label = "False positive rate"
      add_plot.set_xlabel(label, fontsize=12)
      label = "True positive rate"
      add_plot.set_ylabel(label, fontsize=12)
      title = "ROC curve (area = %.4f)"%self._auroc
      add_plot.set_title(title, fontsize=14)
      add_plot.set_aspect(1./add_plot.get_data_ratio())
      add_plot.grid(True)

      print "Assessment"

      add_plot = fig.add_subplot(1, 2, 2)

      gfp = GFP.objects.filter(id=ID)
      gfp = gfp[0]

      path = ""
      path = BASE_DIR + '/data/' + gfp.Type + '/' + gfp.Network + \
	       '/GO_Performance_1_24_2015.mat'
      goperf = loadmat(path, squeeze_me=False, mat_dtype=True)
      goperf = goperf['perf']

      path = ""
      path = BASE_DIR + '/data/' + gfp.Type + '/' + gfp.Network + \
	       '/Phenocarta_Performance_1_24_2015.mat'
      cartaperf = loadmat(path, squeeze_me=False, mat_dtype=True)
      cartaperf = cartaperf['perf']

      path = ""
      path = BASE_DIR + '/data/' + gfp.Type + '/' + gfp.Network + \
	       '/KEGG_Performance_1_24_2015.mat'
      keggperf = loadmat(path, squeeze_me=False, mat_dtype=True)
      keggperf = keggperf['perf']

      add_plot.hist(goperf, bins=linspace(0.4, 1.0, 51), color='forestgreen', alpha=0.5, label="GO")
      add_plot.hist(cartaperf, bins=linspace(0.4, 1.0, 51), color='dodgerblue', alpha=0.5, label="Phenocarta")
      add_plot.hist(keggperf, bins=linspace(0.4, 1.0, 51), color='indianred', alpha=0.5, label="KEGG")
      add_plot.plot((self._auroc, self._auroc), (0, 40), color='indianred', linestyle='--')

      add_plot.legend(loc="upper right", fontsize=6)

      label = "AUROC"
      add_plot.set_xlabel(label, fontsize=12)
      label = "# GO terms"
      add_plot.set_ylabel(label, fontsize=12)
      add_plot.set_aspect(1./add_plot.get_data_ratio())
      add_plot.grid(True)

      filename = BASE_DIR + "/figures/ROC_Curve_" + str(ID) + ".png"
    fig.savefig(filename, format="png", dpi=300)
    del fig

  def plotNodeDegree(self, gfpID):
    (x, y) = self._network.shape
    network = self._network
    network[range(x), range(x)] = 0.0 
    nodedegreenet = network.sum(axis=1)
    indecesann = where(self._annotation == 1)[0]
    (numgenes,) = indecesann.shape
    nodedegree = nodedegreenet[indecesann]

    subnet = network[indecesann,:]
    subnet = subnet[:,indecesann]
    nodedegreesubnet = subnet.sum(axis=1)

    nodedegree = asarray(nodedegree).ravel()
    nodedegreenet = asarray(nodedegreenet).ravel()
    nodedegreesubnet = asarray(nodedegreesubnet).ravel()

    fig = figure()
    fig.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.9, hspace=0.6, wspace=0.4)

    add_plot = fig.add_subplot(1, 2, 1)
    add_plot.set_ylim([0.0, max(nodedegreesubnet)+0.5])

    xx = range(0, int(max(nodedegree)))
    expected = []
    for value in xx:
      expected.append((value * numgenes)/float(x))

    nodedegreescaled = []
    for nd in nodedegree:
      nodedegreescaled.append((nd * numgenes)/float(x))
    add_plot.set_xlim([0.0, max(nodedegreescaled)+1.0])

    identity = arange(0, max(nodedegreesubnet), 0.1)

    add_plot.plot(identity, identity, color='0.5')
    add_plot.scatter(nodedegreescaled, nodedegreesubnet, color='forestgreen', marker='o', alpha=.5, s=10)

    label = "Expected connectivity"
    add_plot.set_xlabel(label, fontsize=12)
    label = "Actual connectivity"
    add_plot.set_ylabel(label, fontsize=12, multialignment='center')

    add_plot.set_aspect(1./add_plot.get_data_ratio())
    add_plot.grid(True)

    add_plot = fig.add_subplot(1, 2, 2)

    residuals = []
    i = 0
    for nd in nodedegreesubnet:
      expected = (nodedegreenet[i] * numgenes) / float(x)
      residuals.append(nd - expected)
      i += 1

    print "residuals: ", residuals
    if sum(residuals) > 0:

      add_plot.set_xlim([min(residuals), max(residuals)])

      pdf = gaussian_kde(residuals)
      x = linspace(min(residuals), max(residuals), 100)
      add_plot.plot(x, pdf(x), 'forestgreen')

    label = "Residual expected connectivity"
    add_plot.set_xlabel(label, fontsize=12)
    label = "Probability amplitude"
    add_plot.set_ylabel(label, fontsize=12, multialignment='center')

    add_plot.set_aspect(1./add_plot.get_data_ratio())
    add_plot.grid(True)

    filename = BASE_DIR + "/figures/Node_Degree_" + str(gfpID) + ".png"
    fig.savefig(filename, format="png", dpi=300)
    del fig

  def plotNetwork(self, gfpID, mapper="Mapper_Symbol_Entrez.pkl"):
    proteinidann = []
    symbolann = []
    annotation = NovelCandidates.objects.filter(gfp=gfpID, Annotation=1, CV=-1)
    for protein in annotation:
      proteinidann.append(protein.ProteinID)
      symbolann.append(protein.Symbol)
    indecesann = where(self._proteinid == proteinidann)[0]

    novel = NovelCandidates.objects.filter(gfp=gfpID, Annotation=0, CV=-1)
    proteinidnovel = []
    symbolnovel = []
    for protein in novel:
      proteinidnovel.append(protein.ProteinID)
      symbolnovel.append(protein.Symbol)
    indecesnovel = where(self._proteinid == proteinidnovel)[0]

    gfp = GFP.objects.filter(id=gfpID)
    gfp = gfp[0]

    if gfp.Type == 'Co':
      network = self._network
      (x, y) = network.shape
      network[range(x), range(x)] = 0.0 
      maxweights = network.max(axis=0)
      threshold = max(maxweights) - ((max(maxweights) * 2)/float(100))
      network[network<threshold] = 0.0
      network[network==threshold] = 0.0
      network[network>threshold] = 1.0

      subnet = network[indecesann,:]
      subnet = subnet[:, indecesnovel]

      check = subnet.sum(axis=1)
      check = asarray(check).ravel()
      indecescheck = where(check > 0.0)[0]
      indecesanncon = indecesann[indecescheck]
      proteinidanncon = asarray(self._proteinid[indecesanncon]).ravel()

      indeces = append(indecesanncon, indecesnovel)

      count = network.sum(axis=0)
      count = sum(count)

      G = nx.from_numpy_matrix(network)

      indecespathscan = {}
      count = 0
      for source in indeces:
        for target in indeces:
          if source != target:
            i = 0
            for path in nx.all_simple_paths(G, source, target, cutoff=2):
              for node in path:
                if node not in indeces:
                  if node in indecespathscan.keys():
                    indecespathscan[node] += 1
                  else:
                    indecespathscan[node] = 1
                if len(indecespathscan.keys()) > 50:
                  break
              i += 1
              if i == 15 or len(indecespathscan.keys()) > 50:
                break
        count += 1
        if count == 1:
          break

      c = Counter(indecespathscan)
      indecespath = {}
      for node, counts in c.most_common(20):
        indecespath[node] = counts	

      indeces = append(indeces, indecespath.keys())
      indeces = indeces.astype(int)

      subnet = network[indeces,:]
      subnet = subnet[:,indeces]

      fig = figure()
      add_plot = fig.add_subplot(1, 1, 1)

      g = nx.from_numpy_matrix(subnet)
      pos = nx.graphviz_layout(g, prog='neato')

      nodelistann = []
      nodelistannsymbol = []
      nodelistnovel = []
      nodelistnovelsymbol = []
      nodelistpaths = []
      nodelistpathssymbol = []

      i = 0
      for index in indeces:
        if index in indecesanncon:
          nodelistann.append(i)
          nodelistannsymbol.append(symbolann[proteinidann.index(self._proteinid[index])])
        if index in indecesnovel:
          nodelistnovel.append(i)
          nodelistnovelsymbol.append(symbolnovel[proteinidnovel.index(self._proteinid[index])])
        if index in indecespath.keys():
          nodelistpaths.append(i)
          symbols = GeneSymbol.objects.filter(ProteinID=int(self._proteinid[index]))
          symbol = ""
          for symbol in symbols:
            symbol = symbol.Symbol
            break
          nodelistpathssymbol.append(symbol)
        i += 1

      nx.draw_networkx_nodes(g, pos, nodelist=nodelistpaths, node_color='grey', alpha=0.2, node_size=10)
      nx.draw_networkx_nodes(g, pos, nodelist=nodelistann, node_color='dodgerblue', alpha=0.65, node_size=40)
      nx.draw_networkx_nodes(g, pos, nodelist=nodelistnovel, node_color='indianred', alpha=0.65, node_size=60)
      nx.draw_networkx_edges(g, pos, alpha=0.2, edge_color='grey')

      edgelist = []
      for edge in g.edges_iter(nodelistann+nodelistnovel):
        if edge[0] in nodelistann + nodelistnovel and edge[1] in nodelistann + nodelistnovel:
          edgelist.append(edge)
				
	
      nx.draw_networkx_edges(g, pos, edgelist=edgelist, edge_color="forestgreen", alpha=0.35)
	
      labels = {}
      for node in g.nodes():
        if node in nodelistann:
          index = nodelistann.index(node)
          labels[node] = nodelistannsymbol[index]
      nx.draw_networkx_labels(g, pos, labels, font_size=6, font_color="black")

      labels = {}
      for node in g.nodes():
        if node in nodelistnovel:
          index = nodelistnovel.index(node)
          labels[node] = nodelistnovelsymbol[index]
      nx.draw_networkx_labels(g, pos, labels, font_size=6, font_color="black")
	
      labels = {}
      for node in g.nodes():
        if node in nodelistpaths:
          index = nodelistpaths.index(node)
          labels[node] = nodelistpathssymbol[index]
      nx.draw_networkx_labels(g, pos, labels, font_size=4, font_color="black")
      axis("off")

      filename = BASE_DIR + "/figures/Network_" + str(gfpID) + ".png"
      fig.savefig(filename, format="png", dpi=350)
      del fig	
    else:
      (x, y) = self._network.shape
      self._network[range(x), range(x)] = 0.0
      subnet = self._network[indecesann,:]
      subnet = subnet[:, indecesnovel]

      check = subnet.sum(axis=1)
      check = asarray(check).ravel()
      indecescheck = where(check > 0.0)[0]
      indecesanncon = indecesann[indecescheck]
      proteinidanncon = asarray(self._proteinid[indecesanncon]).ravel()

      indeces = append(indecesanncon, indecesnovel)

      G = nx.from_numpy_matrix(self._network)

      indecespathscan = {}
      count = 0
      for source in indeces:
        for target in indeces:
          if source != target:
            i = 0
            for path in nx.all_simple_paths(G, source, target, cutoff=2):
              for node in path:
                if node not in indeces:
                  if node in indecespathscan.keys():
                    indecespathscan[node] += 1
                  else:
                    indecespathscan[node] = 1
              i += 1
              if i == 15 or len(indecespathscan.keys()) > 50:
                break
        count += 1
        if count == 1:
          break

      c = Counter(indecespathscan)
      indecespath = {}
      for node, counts in c.most_common(20):
        indecespath[node] = counts

      indeces = append(indeces, indecespath.keys())
      indeces = indeces.astype(int)

      subnet = self._network[indeces,:]
      subnet = subnet[:,indeces]

      fig = figure()
      add_plot = fig.add_subplot(1, 1, 1)

      g = nx.from_numpy_matrix(subnet)
      pos = nx.graphviz_layout(g, prog='neato')

      nodelistann = []
      nodelistannsymbol = []
      nodelistnovel = []
      nodelistnovelsymbol = []
      nodelistpaths = []
      nodelistpathssymbol = []

      i = 0
      for index in indeces:
        if index in indecesanncon:
          nodelistann.append(i)
          nodelistannsymbol.append(symbolann[proteinidann.index(self._proteinid[index])])
        if index in indecesnovel:
          nodelistnovel.append(i)
          nodelistnovelsymbol.append(symbolnovel[proteinidnovel.index(self._proteinid[index])])
        if index in indecespath.keys():
          nodelistpaths.append(i)
          symbols = GeneSymbol.objects.filter(ProteinID=int(self._proteinid[index]))
          symbol = ""
          for symbol in symbols:
            symbol = symbol.Symbol
            break
          nodelistpathssymbol.append(symbol)
			
        i += 1

      nx.draw_networkx_nodes(g, pos, nodelist=nodelistpaths, node_color='grey', alpha=0.2, node_size=10)
      nx.draw_networkx_nodes(g, pos, nodelist=nodelistann, node_color='dodgerblue', alpha=0.65, node_size=40)
      nx.draw_networkx_nodes(g, pos, nodelist=nodelistnovel, node_color='indianred', alpha=0.65, node_size=60)
      nx.draw_networkx_edges(g, pos, alpha=0.2, edge_color='grey')

      edgelist = []
      for edge in g.edges_iter(nodelistann+nodelistnovel):
        if edge[0] in nodelistann + nodelistnovel and edge[1] in nodelistann + nodelistnovel:
          edgelist.append(edge)
      nx.draw_networkx_edges(g, pos, edgelist=edgelist, edge_color="forestgreen", alpha=0.35)

      labels = {}
      for node in g.nodes():
        if node in nodelistann:
          index = nodelistann.index(node)
          labels[node] = nodelistannsymbol[index]
      nx.draw_networkx_labels(g, pos, labels, font_size=6, font_color="black")

      labels = {}
      for node in g.nodes():
        if node in nodelistnovel:
          index = nodelistnovel.index(node)
          labels[node] = nodelistnovelsymbol[index]
      nx.draw_networkx_labels(g, pos, labels, font_size=6, font_color="black")

      labels = {}
      for node in g.nodes():
        if node in nodelistpaths:
          index = nodelistpaths.index(node)
          labels[node] = nodelistpathssymbol[index]
      nx.draw_networkx_labels(g, pos, labels, font_size=4, font_color="black")
      axis("off")

      filename = BASE_DIR + "/figures/Network_" + str(gfpID) + ".png"
      print filename
      fig.savefig(filename, format="png", dpi=350)
      del fig	
	
  @transaction.atomic
  def save(self, gfpID):

    ranks = self._standRank(self._scores)
    i = 0
    predictions = []
    for protein in self._proteinid:
      predictions.append(Prediction(gfp_id=gfpID, CV=-1, ProteinID=protein, \
                                    Score=self._scores[i], \
                                    Rank=ranks[i], \
                                    Annotation=self._annotation[i]))
      i += 1
    try:
      Prediction.objects.bulk_create(predictions, batch_size=500)
    except Exception as e:
      print "Exception= ", e
      print type(e)
      print "protein= ", protein
      print e.args

    ranks = self._standRank(self._foldscores)
    i = 0
    predictions = []
    for protein in self._foldproteins:
      predictions.append(Prediction(gfp_id=gfpID, ProteinID=protein, \
                                    Score=self._foldscores[i], \
                                    Rank=ranks[i], \
                                    Annotation=self._foldlabels[i], CV=3))
      i += 1
    try:
      Prediction.objects.bulk_create(predictions, batch_size=500)
    except Exception as e:
      log.debug("RUNNING: save bulk_create exception %s", str(e))
      log.debug("RUNNING: save bulk_create exception %s", str(type(e)))
      #print e.args

    per = Performance()
    per.save(gfpID, self._auroc)
    del per

    predictions = Prediction.objects.filter(gfp=gfpID, Annotation=0, CV=-1).order_by('Score')[0:20]

    j = 1
    for prediction in predictions:
      symbols = GeneSymbol.objects.filter(ProteinID=prediction.ProteinID)
      symbol = ""

      for symbol in symbols:
        symbol = symbol.Symbol
        break
      NovelCandidates.objects.get_or_create(gfp_id=gfpID, Symbol=symbol, \
                                            ProteinID=prediction.ProteinID, \
                                            Score=prediction.Score, \
                                            Rank=prediction.Rank, CV=-1, \
                                            Annotation=0)

    predictions = Prediction.objects.filter(gfp=gfpID, Annotation=0, CV=3).order_by('Score')[0:20]

    j = 1
    for prediction in predictions:
      symbols = GeneSymbol.objects.filter(ProteinID=prediction.ProteinID)
      symbol = ""

      for symbol in symbols:
        symbol = symbol.Symbol
        break
      NovelCandidates.objects.get_or_create(gfp_id=gfpID, Symbol=symbol, \
                                            ProteinID=prediction.ProteinID, \
                                            Score=prediction.Score, \
                                            Rank=prediction.Rank, CV=3, \
                                            Annotation=0)

    log.debug("RUNNING: save plots")
    self.plotROCCurve(gfpID)
    self.plotNetwork(gfpID)
    self.plotNodeDegree(gfpID)

  def enrichment(self, experimentID, gfpID):

    gfp = GFP.objects.filter(id=gfpID)
    gfp = gfp[0]

    Type = ""
    if gfp.Network == "Co":
      print "GFP network: ", gfp.Network
      network = CoAnnotation.objects.filter(Geoid=gfp.Network)
      print "len network: ", len(network)
      print network
      network = network[0]
      Type = network.Type
      print "Type: ", Type
    else:
      Type = gfp.Type
      print "Type: ", Type 

    # GO enrichment
    enrich = EnrichmentAnalysis(experimentID, gfpID)
    enrich.loadProfile("GOID_Human_1_24_2015.mat", \
                       "Gene2GO_Human_1_24_2015.mat")

    if not Enrichment.objects.filter(experiment_id=experimentID, Profile="GO").exists():
      enrich.run(profile="GO")
    if Type == 'sc RNA Seq':
      enrich.runGFPscRNASeq(profile="GO")
    else:
      enrich.runGFP(profile="GO")
    del enrich
    
  def aggregate(self, gfpIDs, experimentID):
    proteinid = Prediction.objects.filter(gfp=gfpIDs[0], CV=-1).order_by("ProteinID").values_list('ProteinID', flat=True)

    label = Prediction.objects.filter(gfp=gfpIDs[0], CV=-1).order_by("ProteinID").values_list('Annotation', flat=True)
    label = asarray(label)

    y = len(Prediction.objects.filter(gfp=gfpIDs[0], CV=-1).order_by("ProteinID").values_list('Score', flat=True))
    prediction = zeros((len(gfpIDs), y))

    proteinidcv = Prediction.objects.filter(gfp=gfpIDs[0], CV=3).order_by("ProteinID").values_list('ProteinID', flat=True)

    labelcv = Prediction.objects.filter(gfp=gfpIDs[0], CV=3).order_by("ProteinID").values_list('Annotation', flat=True)
    labelcv = asarray(labelcv)

    y = Prediction.objects.filter(gfp=gfpIDs[0], CV=3).order_by("ProteinID").count()
    performance = zeros((len(gfpIDs), y))

    i = 0
    for gfpid in gfpIDs:
      predictor = Prediction.objects.filter(gfp=gfpid, CV=-1).order_by("ProteinID").values_list('Score', flat=True)
      predictor = asarray(predictor).ravel()
      prediction[i, :] = asarray(predictor)
	
      performer = Prediction.objects.filter(gfp=gfpid, CV=3).order_by("ProteinID").values_list('Score', flat=True)
      performer = asarray(performer).ravel()
      performance[i, :] = asarray(performer)

      i += 1

    aggpred = prediction.mean(axis=0)
    aggpred = asarray(aggpred.ravel())

    rankpred = self._standRank(aggpred)
	
    aggperf = performance.mean(axis=0)
    aggperf = asarray(aggperf.ravel())

    rankperf = self._standRank(aggperf)

    i = 0
    flag = AggregatePrediction.objects.filter(experiment=experimentID).exists()
    if flag:
      predictions = AggregatePrediction.objects.filter(experiment=experimentID, CV=-1).order_by("ProteinID")
      i = 0
      for prediction in predictions:
        prediction.Score = aggpred[i]
        prediction.Rank = rankpred[i]
        i += 1

      try:
        bulk_update(predictions, batch_size=500)
      except Exception as e:
        log.debug("RUNNING: save bulk_update exception %s", str(e))
        log.debug("RUNNING: save bulk_update exception %s", str(type(e)))
    else:
      predictions = []
      i = 0
      for protein in proteinid:
        predictions.append(AggregatePrediction(experiment_id=experimentID, CV=-1, \
                                               ProteinID=protein, Annotation=label[i], \
                                               Score=aggpred[i], Rank=rankpred[i]))
        i += 1
      try:
        AggregatePrediction.objects.bulk_create(predictions, batch_size=500)
      except Exception as e:
        log.debug("RUNNING: save bulk_create exception %s", str(e))
        log.debug("RUNNING: save bulk_create exception %s", str(type(e)))

    if flag:
      predictions = AggregatePrediction.objects.filter(experiment=experimentID, CV=3).order_by("ProteinID")
      i = 0
      for prediction in predictions:
        prediction.Score = aggperf[i]
        prediction.Rank = rankperf[i]
        i += 1
      try:
        bulk_update(predictions, batch_size=500)
      except Exception as e:
        log.debug("RUNNING: save bulk_update exception %s", str(e))
        log.debug("RUNNING: save bulk_update exception %s", str(type(e)))
    else:
      predictions = []
      i = 0
      for protein in proteinidcv:
        predictions.append(AggregatePrediction(experiment_id=experimentID, CV=3, \
                                               ProteinID=protein, Annotation=labelcv[i], \
                                               Score=aggperf[i], Rank=rankperf[i]))
        i += 1
      try:
        AggregatePrediction.objects.bulk_create(predictions, batch_size=500)
      except Exception as e:
        log.debug("RUNNING: save bulk_create exception %s", str(e))
        log.debug("RUNNING: save bulk_create exception %s", str(type(e)))

    self._auroc = self.AUROC(labelcv, aggperf)

    # check is exists
    if AggregatePerformance.objects.filter(experiment=experimentID).exists():
      AggregatePerformance.objects.filter(experiment=experimentID).update(AUROC=self._auroc)
    else:
      per = AggregatePerformance()
      per.save(experimentID, self._auroc)
      del per

    # representation bias
    nodedegree = zeros(self._proteinid.shape)
    nodedegree = nodedegree.ravel()
    print "ND: ", nodedegree.shape
    
    for gfp in gfpIDs:
      path = ""
      if gfp.Type == 'PPI':
        path = BASE_DIR + "/data/"+ gfp.Type + "/" + gfp.Network + "/ppiNetwork_Human_1_24_2015.mat"
      elif gfp.Type == 'Co':
        path = BASE_DIR + "/data/"+ gfp.Type + "/" + gfp.Network + "/coNetwork_Human_1_24_2015.mat" 
      elif gfp.Type == 'SM':
        path = BASE_DIR + "/data/"+ gfp.Type + "/" + gfp.Network + "/smNetwork_Human_1_24_2015.mat"
      else:
        log.debug("Invalid network - type")

      network = loadmat(path, squeeze_me=False, mat_dtype=True)
      network = network['network']
      if issparse(network):
        network = network.todense()

      nodedegreenet = network.sum(axis=1)
      nodedegreenet = asarray(nodedegreenet).ravel()

      nodedegree += nodedegreenet

    nodedegree = rankdata(nodedegree)
    (x, y) = self._proteinid.shape
    nodedegree /= float(x)

    print "ND: ", nodedegree.shape

    rrep = []

    path = BASE_DIR + "/data/representation_PPI.mat"
    predictor = loadmat(path, squeeze_me=False, mat_dtype=True)
    representation_PPI = predictor['predictor']
    representation_PPI = asarray(representation_PPI).ravel()
    print spearmanr(representation_PPI, nodedegree)
    (ppir, ppip) = spearmanr(representation_PPI, nodedegree)
    rrep.append(ppir)

    path = BASE_DIR + "/data/representation_SM.mat"
    predictor = loadmat(path, squeeze_me=False, mat_dtype=True)
    representation_SM = predictor['predictor']
    representation_SM = asarray(representation_SM).ravel()
    print spearmanr(representation_SM, nodedegree)
    (smr, smp) = spearmanr(representation_SM, nodedegree)
    rrep.append(smr)

    path = BASE_DIR + "/data/representation_scRNASeq.mat"
    predictor = loadmat(path, squeeze_me=False, mat_dtype=True)
    representation_scRNASeq = predictor['predictor']
    representation_scRNASeq = asarray(representation_scRNASeq).ravel()
    print spearmanr(representation_scRNASeq, nodedegree)
    (scRNASeqr, scRNASeqp) = spearmanr(representation_scRNASeq, nodedegree)
    rrep.append(scRNASeqr)

    path = BASE_DIR + "/data/representation_RNASeq.mat"
    predictor = loadmat(path, squeeze_me=False, mat_dtype=True)
    representation_RNASeq = predictor['predictor']
    representation_RNASeq = asarray(representation_RNASeq).ravel()
    print spearmanr(representation_RNASeq, nodedegree)
    (RNASeqr, RNASeqp) = spearmanr(representation_RNASeq, nodedegree)
    rrep.append(RNASeqr)

    path = BASE_DIR + "/data/representation_Micro1.mat"
    predictor = loadmat(path, squeeze_me=False, mat_dtype=True)
    representation_micro = predictor['predictor']
    representation_micro = asarray(representation_micro).ravel()
    print spearmanr(representation_micro, nodedegree)
    (micror, microp) = spearmanr(representation_micro, nodedegree)
    rrep.append(micror)

    self._TPR_FPR(labelcv, aggperf)
    self.plotROCCurve(experimentID, aggregate=True, rep=rrep)

  # Factory method
  @staticmethod
  def getAlgorithm(name):
    if name == "NV":
      return NV()
    elif name == "LR":
      return LR()
    elif name == "RWR":
      return RWR()
    elif name == "SGD":
      return SGD()
    elif name == "PA":
      return PA()
    elif name == "RP":
      return RP()
    else:
      return None


class NV(Algorithm):
  def __init__(self):
    Algorithm.__init__(self)
    self.__name = "NV"
  def loadData(self, gfpID):
    Algorithm.loadData(self, gfpID)
  def _fit(self, annotation):
    sumin = dot(self._network, annotation)
    sumall = tile(sum(self._network), (1, 1))
    sumall = sumall.transpose()
    scores = sumin/sumall
    del sumin
    del sumall

    scores = subtract(ones_like(scores), scores)
    return scores

  def run(self, nFold=3):
    log.debug("NV: run")
    (numx, numy) = self._network.shape

    pp = permutation(numx)
    self._scores = self._fit(self._annotation)
    auroc = self.AUROC(self._annotation, self._scores)

    fold = 0
    offset = 0
    meanroc = []
    while fold < nFold:
      log.debug("NV: ___ fold= %d ___" % fold)
      lastelem = int(min(numx, offset+floor(numx/nFold)))

      ix = []
      for index in pp[offset+1:lastelem]:
        ix.append(index)

      offset = lastelem

      labeltmp = []
      for value in self._annotation:
        labeltmp.append(float(value))

      for index in ix:
        labeltmp[index] = 0

      labeltmp = asarray(labeltmp)
      labeltmp = labeltmp.reshape(labeltmp.shape[0], 1)

      try:
        scores = self._fit(labeltmp)
      except Exception as e:
        log.debug("NV: fit exception %s", str(e))
        log.debug("NV: fit exception %s", str(type(e)))
        break

      score = []
      label = []
      protein = []
      for index in ix:
        score.append(float(scores[index]))
        label.append(self._annotation[index])
        protein.append(self._proteinid[index])

        self._foldlabels.append(int(self._annotation[index]))
        self._foldscores.append(float(scores[index]))
        self._foldproteins.append(int(self._proteinid[index]))

      auroc = self.AUROC(label, score)
      log.debug("AUROC= %.4f" % auroc)

      meanroc.append(auroc)

      fold += 1

    self._auroc = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
    self._TPR_FPR(self._foldlabels, self._foldscores)


class LR(Algorithm):
  def __init__(self):
    Algorithm.__init__(self)
    self.__name = "LR"
  def loadData(self, gfpID):
    Algorithm.loadData(self, gfpID)
  def run(self, nFold=3, dual=True, penalty='l2', class_weight='auto', C=1.0, tol=0.0001):
    log.debug("LR: run")
    (numx, numy) = self._network.shape

    pp = permutation(numx)
		
    model = LogisticRegression(penalty=penalty, dual=dual, class_weight=class_weight, C=C, tol=tol)
    model.fit(self._network, self._annotation.ravel())
    scores = model.predict_proba(self._network)

    self._scores = scores[:,0]

    fold = 0
    offset = 0
    meanroc = []
    labelIx = range(numx)
    while fold < nFold:
      log.debug("NV: ___ fold= %d ___" % fold)
      lastelem = int(min(numx, offset+floor(numx/nFold)))

      ix = []
      for index in pp[offset+1:lastelem]:
        ix.append(index)

      offset = lastelem

      labeltmp = []
      for value in self._annotation:
        labeltmp.append(float(value))

      for index in ix:
        labeltmp[index] = 0

      model = LogisticRegression(penalty=penalty, dual=dual, class_weight=class_weight, C=C, tol=tol)
      model.fit(self._network, labeltmp)
      scores = model.predict_proba(self._network)
      scores = scores[:,0]

      score = []
      label = []
      protein = []
      for index in ix:
        score.append(float(scores[index]))
        label.append(self._annotation[index])
        protein.append(self._proteinid[index])

        self._foldlabels.append(int(self._annotation[index]))
        self._foldscores.append(float(scores[index]))
        self._foldproteins.append(int(self._proteinid[index]))

      auroc = self.AUROC(label, score)
      log.debug("AUROC= %.4f" % auroc)

      meanroc.append(auroc)

      fold += 1

    self._auroc = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
    self._TPR_FPR(self._foldlabels, self._foldscores)


class RWR(Algorithm):
  def __init__(self):
    Algorithm.__init__(self)
    self.__name = "RWR"
    self.__epsilon = 0.0
    self.__theta = 0.0
  def _model(self, epsilon=0.01, theta=0.8):
    self.__epsilon = epsilon
    self.__theta = theta
  def _fit(self, network, annotation):
    # Initialisation
    (x, y) = network.shape
    p_t1 = zeros(x)
    Vm = len([x for x in annotation if x == 1])
    indeces = where(annotation == 1)[0]
    p_t1[indeces] = 1/float(Vm)
    p_t0 = copy(p_t1)

    i = 0
    d1 = 0.0
    while 1:
      p_t2 = (1  - self.__theta) * dot(self._transition.T, p_t1) + self.__theta * p_t0
      d2 = sum(p_t2-p_t1)

      i += 1
      if i == 100:
        #print "steps"
        return (1 - p_t2)
      elif d2 == 0.0:
        return (1 - p_t2)
      elif (d2+d1) == 0.0:
        return (1 - p_t2)
      d1 = sum(p_t2-p_t1)	
      p_t1 = p_t2
  def loadData(self, gfpID):
    Algorithm.loadData(self, gfpID)
  def run(self, nFold=3):
    log.debug("RWR: run")
    (numx, numy) = self._network.shape

    pp = permutation(numx)
		
    self._model()
		
    self._scores = self._fit(self._network, self._annotation.ravel())

    fold = 0
    offset = 0
    meanroc = []
    labelIx = range(numx)
    while fold < nFold:
      log.debug("RWR: ___ fold= %d ___" % fold)
      lastelem = int(min(numx, offset+floor(numx/nFold)))

      ix = []
      for index in pp[offset+1:lastelem]:
        ix.append(index)

      offset = lastelem

      labeltmp = []
      for value in self._annotation:
        labeltmp.append(float(value))

      for index in ix:
        labeltmp[index] = 0

      labeltmp = asarray(labeltmp)
      labeltmp = labeltmp.ravel()

      self._model()
      scores = self._fit(self._network, labeltmp)

      score = []
      label = []
      protein = []
      for index in ix:
        score.append(float(scores[index]))
        label.append(int(self._annotation[index]))
        protein.append(int(self._proteinid[index]))

        self._foldlabels.append(int(self._annotation[index]))
        self._foldscores.append(float(scores[index]))
        self._foldproteins.append(int(self._proteinid[index]))

      auroc = self.AUROC(label, score)
      log.debug("AUROC= %.4f" % auroc)

      meanroc.append(auroc)

      fold += 1

    self._auroc = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
    auroc = self.AUROC(self._foldlabels, self._foldscores)
    self._TPR_FPR(self._foldlabels, self._foldscores)


class SGD(Algorithm):
  def __init__(self):
    Algorithm.__init__(self)
    self.__name = "SGD"
  def loadData(self, gfpID):
    Algorithm.loadData(self, gfpID)
  def run(self, nFold=3, iter=10, verbose=1, loss='modified_huber', penalty='l2', class_weight='auto', \
          shuffle=True):
    log.debug("SGD: run")
    (numx, numy) = self._network.shape

    pp = permutation(numx)
		
    print "0s=", len([x for x in self._annotation if x == 0])
    print "1s=", len([x for x in self._annotation if x == 1])

    model = SGDClassifier(loss=loss, class_weight=class_weight, penalty=penalty, n_iter=iter, shuffle=shuffle, \
								verbose=verbose)
    model.fit(self._network, self._annotation.ravel())
    scores = model.predict_proba(self._network)

    self._scores = scores[:,0]

    fold = 0
    offset = 0
    meanroc = []
    labelIx = range(numx)
    while fold < nFold:
      log.debug("SGD: ___ fold= %d ___" % fold)
      lastelem = int(min(numx, offset+floor(numx/nFold)))

      ix = []
      for index in pp[offset+1:lastelem]:
        ix.append(index)

      offset = lastelem

      labeltmp = []
      for value in self._annotation:
        labeltmp.append(float(value))

      for index in ix:
        labeltmp[index] = 0

      model = SGDClassifier(loss=loss, class_weight=class_weight, penalty=penalty, n_iter=iter, \
                            shuffle=shuffle, verbose=verbose)
      model.fit(self._network, labeltmp)
      scores = model.predict_proba(self._network)
      scores = scores[:,0]

      score = []
      label = []
      protein = []
      for index in ix:
        score.append(float(scores[index]))
        label.append(self._annotation[index])
        protein.append(self._proteinid[index])

        self._foldlabels.append(int(self._annotation[index]))
        self._foldscores.append(float(scores[index]))
        self._foldproteins.append(int(self._proteinid[index]))
			
      auroc = self.AUROC(label, score)
      log.debug("AUROC= %.4f" % auroc)

      meanroc.append(auroc)
      fold += 1

    self._auroc = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
    auroc = self.AUROC(self._foldlabels, self._foldscores)

    self._TPR_FPR(self._foldlabels, self._foldscores)


class PA(Algorithm):
  def __init__(self):
    Algorithm.__init__(self)
    self.__name = "PA"
  def loadData(self, gfpID):
    Algorithm.loadData(self, gfpID)
  def _convertScore(self, scores):
    tmp = []
    for score in scores:
      tmp.append(fabs(score))
    mag = max(tmp)
    del tmp

    scorereturn = []
    for score in scores:
      value = score/float(mag)
      scorereturn.append(1 - ((value/float(2))+0.5))
    return scorereturn
  def run(self, nFold=3, loss='hinge', iter=10, verbose=1):
    log.debug("PA: run")
    (numx, numy) = self._network.shape

    pp = permutation(numx)
		
    model = PassiveAggressiveClassifier(loss=loss, n_iter=iter, verbose=verbose)
    model.fit(self._network, self._annotation.ravel())
    scores = model.decision_function(self._network)
    self._scores = self._convertScore(scores)

    fold = 0
    offset = 0
    meanroc = []
    labelIx = range(numx)
    while fold < nFold:
      log.debug("NV: ___ fold= %d ___" % fold)
      lastelem = int(min(numx, offset+floor(numx/nFold)))

      ix = []
      for index in pp[offset+1:lastelem]:
        ix.append(index)

      print lastelem

      offset = lastelem

      labeltmp = []
      for value in self._annotation:
        labeltmp.append(float(value))

      for index in ix:
        labeltmp[index] = 0

      model = PassiveAggressiveClassifier(loss=loss, n_iter=iter, verbose=verbose)
      model.fit(self._network, labeltmp)
      scores = model.decision_function(self._network)
      scores = self._convertScore(scores)

      score = []
      label = []
      protein = []
      for index in ix:
        score.append(float(scores[index]))
        label.append(int(self._annotation[index]))
        protein.append(int(self._proteinid[index]))

        self._foldlabels.append(int(self._annotation[index]))
        self._foldscores.append(float(scores[index]))
        self._foldproteins.append(int(self._proteinid[index]))

      auroc = self.AUROC(label, score)
      log.debug("AUROC= %.4f" % auroc)

      meanroc.append(auroc)

      fold += 1

    self._auroc = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
    auroc = self.AUROC(self._foldlabels, self._foldscores)

    self._TPR_FPR(self._foldlabels, self._foldscores)


class RP(Algorithm):
  def __init__(self):
    Algorithm.__init__(self)
    self.__name = "RP"
    self.__alpha = 0.0

  def _model(self, alpha=0.95):
    self.__alpha = alpha

  def _reweightNet(self):
    (x, y) = self._network.shape

    net = self._network
    net[range(x), range(x)] = 0
    net = (net.T/(net.sum(axis=1) + spacing(0))).T
    return net

  def _fit(self, network, annotation):
    (x, y) = network.shape
    annotation = asarray(annotation).ravel()
    label = asarray(annotation).ravel()
    y_0 = annotation/sum(annotation)

    y_t = asarray(annotation).ravel()
    y_t = y_t/sum(y_t)

    count = 0
    while count < 20:
      help = y_t * self.__alpha
      help = asarray(help.T).ravel()
      y_t1 = y_0 + dot(network, help)
      y_t1 = asarray(y_t1.T).ravel()
      y_t = y_t1
      count += 1

    return 1 - y_t
    
  def loadData(self, gfpID):
    Algorithm.loadData(self, gfpID)

  def run(self, nFold=3):
    log.debug("RP: run")
    (numx, numy) = self._network.shape

    net = self._reweightNet() 

    pp = permutation(numx)
		
    self._model()
    self._scores = self._fit(net, self._annotation)

    fold = 0
    offset = 0
    meanroc = []
    labelIx = range(numx)
    while fold < nFold:
      log.debug("RP: ___ fold= %d ___" % fold)
      lastelem = int(min(numx, offset+floor(numx/nFold)))

      ix = []
      for index in pp[offset+1:lastelem]:
        ix.append(index)

      offset = lastelem

      labeltmp = []
      for value in self._annotation:
        labeltmp.append(float(value))

      for index in ix:
        labeltmp[index] = 0

      labeltmp = asarray(labeltmp)
      labeltmp = labeltmp.ravel()

      self._model()
      scores = self._fit(net, labeltmp)

      score = []
      label = []
      protein = []
      for index in ix:
        score.append(float(scores[index]))
        label.append(int(self._annotation[index]))
        protein.append(int(self._proteinid[index]))

        self._foldlabels.append(int(self._annotation[index]))
        self._foldscores.append(float(scores[index]))
        self._foldproteins.append(int(self._proteinid[index]))

      auroc = self.AUROC(label, score)
      log.debug("AUROC= %.4f" % auroc)

      meanroc.append(auroc)

      fold += 1

    self._auroc = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
    auroc = self.AUROC(self._foldlabels, self._foldscores)
    self._TPR_FPR(self._foldlabels, self._foldscores)


class TestAlgorithm(TestCase):
  def setUp(self):
    self.algo = Algorithm()
    #self.algo.test("algor construct")
  def tearDown(self):
    print "tearDown"
    del self.algo
  def testALR(self):
    algo = Algorithm.getAlgorithm("LR")
    algo.loadData(132)
    #algo.test("check for LR")

    algo.plotNetwork(132)


def suite():
  suite = makeSuite(TestAlgorithm, 'test')
  return suite

if __name__ == '__main__':
  main(defaultTest='suite', argv=['-d'])
