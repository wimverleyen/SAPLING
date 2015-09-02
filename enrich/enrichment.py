
from unittest import TestCase, makeSuite, main

import pickle
import logging
from numpy import zeros, where, asarray, array, empty, reshape
from scipy.stats import hypergeom, mannwhitneyu
from scipy.misc import comb
from scipy.sparse import issparse
from scipy.io import loadmat, savemat
from statsmodels.stats.multitest import multipletests

from saplingdev.settings import BASE_DIR
from GFPdevs.models import GFP, NovelCandidates, Enrichment, EnrichmentGFP, Prediction, CoAnnotation

from bulk_update.helper import bulk_update

log = logging.getLogger(__name__)

class EnrichmentAnalysis:
  def __init__(self, experimentID, gfpID, ProteinID="entrezGeneID_Human_1_24_2015.mat", \
			mapper="Mapper_Symbol_Entrez.pkl"):

    self._experimentid = experimentID
    self._gfpid = gfpID

    self._profile = None
    self._profileid = None

    # ProteinID
    path = BASE_DIR + "/data/" + ProteinID
    proteinid = loadmat(path, squeeze_me=False, mat_dtype=True)
    self._proteinid = proteinid['proteinid']

    # Annotation
    proteinidann = []
    annotation = NovelCandidates.objects.filter(gfp=gfpID, Annotation=1)
    for protein in annotation:
      proteinidann.append(protein.ProteinID)

    indecesann = where(self._proteinid == proteinidann)[0]
    self._annotation = zeros(self._proteinid.shape)
    self._annotation[indecesann] = 1
    self._annotation = asarray(self._annotation)
    self._annotation = self._annotation.ravel()

    print self._annotation

    # Top 20 novel annotations
    proteinidnovel = []
    annotation = NovelCandidates.objects.filter(gfp=gfpID, Annotation=0)
    for protein in annotation:
      proteinidnovel.append(protein.ProteinID)
    indecesnovel = where(self._proteinid == proteinidnovel)[0]

    self._novel = zeros(self._proteinid.shape)
    self._novel[indecesann] = 1
    self._novel[indecesnovel] = 1
    self._novel = asarray(self._novel).ravel()
    
    proteinid = Prediction.objects.filter(gfp=gfpID, CV=-1).order_by("ProteinID").values_list('ProteinID', \
													flat=True)
    self._predictor = Prediction.objects.filter(gfp=gfpID, CV=-1).order_by("ProteinID").values_list('Score', \
													flat=True)
    self._predictor = asarray(self._predictor).ravel()

  def __del__(self):
    del self._experimentid
    del self._gfpid
    del self._proteinid
    del self._annotation
    del self._profileid
    del self._profile

  def prepareProfile(self, profileid, profile):
    path = BASE_DIR + "/data/" + profileid
    profileid = loadmat(path, squeeze_me=False, mat_dtype=True)
    profileid = profileid['GOID']
    savemat(path, {'profileid':profileid})

    path = BASE_DIR + "/data/" + profile
    profile = loadmat(path, squeeze_me=False, mat_dtype=True)
    profile = profile['go2geneMatrix']
    savemat(path, {'profile':profile})

  def loadProfile(self, profileid, profile):
    path = BASE_DIR + "/data/" + profileid
    profileid = loadmat(path, squeeze_me=False, mat_dtype=True)
    self._profileid = profileid['profileid']

    path = BASE_DIR + "/data/" + profile
    profile = loadmat(path, squeeze_me=False, mat_dtype=True)
    self._profile = profile['profile']
    if issparse(self._profile):
      self._profile = self._profile.todense()
  
  """  
  def runscRNASeq(self, profile="GO"):
    log.debug("Run sc RNA Seq: %s profile" % profile)

    # Load network
    gfp = GFP.objects.filter(id=self._gfpid)
    log.debug("Run sc RNA Seq: network %s" % gfp[0].Network)
    print "Run sc RNA Seq: network %s" % gfp[0].Network

    path = BASE_DIR + "/data/Co/" + gfp[0].Network + "/coNetwork_Human_1_24_2015.mat" 
    network = loadmat(path, squeeze_me=False, mat_dtype=True)
    network = network['network']
    (x, y) = network.shape
    if issparse(network):
      network = network.todense()

    # check for genes with connection in the network
    network[range(x), range(x)] = 0
    nodedegree = network.sum(axis=0)
    nodedegree = asarray(nodedegree).ravel()
    indeces = where(nodedegree > 0)[0]
    indeceszero = where(nodedegree == 0)[0]
    (N,) = indeces.shape
    print "universe: ", N

    (n,) = where(self._annotation == 1)[0].shape
    print "length annotations: ", n

    # Filter
    filter = self._profile[indeces, :]
    filter = filter.sum(axis=0)
    indeces = where(filter > 20)[0]

    prof = self._profile
    prof[indeceszero,:] = 0
    prof = prof[:,indeces]
    profid = self._profileid[indeces]
    
    i = 0
    pvalues= []
    enrich = []
    for id in profid:

      check = where(prof[:, i] == 1)[0]
      check = asarray(check).ravel()
      (K,) = check.shape
      print "K: ", K

      function = asarray(profile[:, i])
      function = function.ravel()
      overlap = function + self._annotation
      (k,) = where(overlap == 2)[0].shape
      print "k: ", k

      if k > 0:
        #print "p= ", self.hypergeometric(N, K, n, k)
        pvalues.append(self.hypergeometric(N, K, n, k))
	enrich.append(Enrichment(experiment_id=self._experimentid, Profile=profile, \
			ProfileID=int(profid[i]), Pvalue= self.hypergeometric(N, K, n, k)))
      else:
       	enrich.append(Enrichment(experiment_id=self._experimentid, Profile=profile, \
			ProfileID=int(profid[i]), Pvalue= 1.0))
      i += 1

    if Enrichment.objects.filter(experiment_id=self._experimentid, Profile=profile).exists():
      enrichment = Enrichment.objects.filter(experiment_id=self._experimentid, Profile=profile).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)
    else:
      Enrichment.objects.bulk_create(enrich, batch_size=500)

      enrichment = Enrichment.objects.filter(experiment_id=self._experimentid, Profile=profile).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)

      i += 1
  """  

  def run(self, profile="GO"):
    log.debug("Run: %s profile" % profile)

    i = 0

    (N, z) = self._proteinid.shape
    (n,) = where(self._annotation == 1)[0].shape

    pvalues= []
    enrich = []
    for id in self._profileid:
      check = where(self._profile[:, i] == 1)[0]
      check = asarray(check).ravel()

      (K,) = check.shape

      function = asarray(self._profile[:, i])
      function = function.ravel()
      overlap = function + self._annotation

      (k,) = where(overlap == 2)[0].shape

      if k > 0:
        pvalues.append(self.hypergeometric(N, K, n, k))
	enrich.append(Enrichment(experiment_id=self._experimentid, Profile=profile, \
			ProfileID=int(self._profileid[i]), Pvalue= self.hypergeometric(N, K, n, k)))
      else:
       	enrich.append(Enrichment(experiment_id=self._experimentid, Profile=profile, \
			ProfileID=int(self._profileid[i]), Pvalue= 1.0))
      i += 1
    if Enrichment.objects.filter(experiment_id=self._experimentid, Profile=profile).exists():
      enrichment = Enrichment.objects.filter(experiment_id=self._experimentid, Profile=profile).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)
    else:
      Enrichment.objects.bulk_create(enrich, batch_size=500)

      enrichment = Enrichment.objects.filter(experiment_id=self._experimentid, Profile=profile).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)

  def runGFPscRNASeq(self, profile="GO"):
    log.debug("RunGFP: %s profile" % profile)
    #print "RunGFP: %s profile" % profile

    # get GFP predictions for the gene list
    gfp = GFP.objects.filter(id=self._gfpid)
    log.debug("Run sc RNA Seq: network %s" % gfp[0].Network)
    #print "Run sc RNA Seq: network %s" % gfp[0].Network

    path = BASE_DIR + "/data/Co/" + gfp[0].Network + "/coNetwork_Human_1_24_2015.mat" 
    network = loadmat(path, squeeze_me=False, mat_dtype=True)
    network = network['network']
    (x, y) = network.shape
    if issparse(network):
      network = network.todense()

    # check for genes with connection in the network
    network[range(x), range(x)] = 0
    nodedegree = network.sum(axis=0)
    nodedegree = asarray(nodedegree).ravel()
    indeces = where(nodedegree > 0)[0]
    indeceszero = where(nodedegree == 0)[0]
    (N,) = indeces.shape
    #print "universe: ", N

    i = 0
    (n,) = where(self._novel == 1)[0].shape
    #print "length annotations: ", n

    # Filter
    filter = self._profile[indeces, :]
    filter = filter.sum(axis=0)
    indeces = where(filter > 20)[0]

    prof = self._profile
    prof[indeceszero,:] = 0
    prof = prof[:,indeces]
    profid = self._profileid[indeces]

    pvalues= []
    enrich = []
    for id in profid:
      check = where(prof[:, i] == 1)[0]
      check = asarray(check).ravel()
      #print check.shape
      (K,) = check.shape
      function = asarray(prof[:, i]).ravel()
      #function = function.ravel()
      #print "novel: ", self._novel.shape
      #print "novel: ", self._novel
      #print "function: ", function.shape
      #print "function: ", function

      overlap = function + self._novel

      (k,) = where(overlap == 2)[0].shape

      if k > 0:
        pvalues.append(self.hypergeometric(N, K, n, k))
	enrich.append(EnrichmentGFP(gfp_id=self._gfpid, Profile=profile, \
			ProfileID=int(profid[i]), Test="hypergeometric", \
			Pvalue= self.hypergeometric(N, K, n, k)))
      else:
       	enrich.append(EnrichmentGFP(gfp_id=self._gfpid, Profile=profile, \
			ProfileID=int(profid[i]), Test="hypergeometric", Pvalue= 1.0))
      i += 1

    if EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile, Test="hypergeometric").exists():
      enrichment = EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)
    else:
      EnrichmentGFP.objects.bulk_create(enrich, batch_size=500)
      enrichment = EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile, Test="hypergeometric").order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)

  def runGFP(self, profile="GO"):
    log.debug("RunGFP: %s profile" % profile)
    print "RunGFP: %s profile" % profile
    # get GFP predictions for the gene list
    score = []
    predictions = Prediction.objects.filter(gfp_id=self._gfpid, CV=-1)
    for prediction in predictions:
      score.append(prediction.Score)

    score = asarray(score)

    i = 0
    pvalues= []
    enrich = []
    for id in self._profileid:
      indeces0 = where(self._profile[:, i] == 0)[0]
      indeces1 = where(self._profile[:, i] == 1)[0]

      x = score[indeces0]
      x = asarray(x).ravel()
      y = score[indeces1]
      y = asarray(y).ravel()

      enrich.append(EnrichmentGFP(gfp_id=self._gfpid, Profile=profile, \
                  ProfileID=int(self._profileid[i]), Test="Mann-Whitney", Pvalue= mannwhitneyu(x, y)[1]))

      i += 1

    if EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile, Test="Mann-Whitney").exists():
      enrichment = EnrichmentGFP.objects.filter(gfp_id=self._gfpid, profile=pro).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)

      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)
    else:
      EnrichmentGFP.objects.bulk_create(enrich, batch_size=500)

      enrichment = EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)

      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)

    i = 0
    (N, z) = self._proteinid.shape
    (n,) = where(self._novel == 1)[0].shape

    pvalues= []
    enrich = []
    for id in self._profileid:
      check = where(self._profile[:, i] == 1)[0]
      check = asarray(check).ravel()
      print check.shape
      (K,) = check.shape
      function = asarray(self._profile[:, i])
      function = function.ravel()
      overlap = function + self._novel

      (k,) = where(overlap == 2)[0].shape

      if k > 0:
        pvalues.append(self.hypergeometric(N, K, n, k))
	enrich.append(EnrichmentGFP(gfp_id=self._gfpid, Profile=profile, \
			ProfileID=int(self._profileid[i]), Test="hypergeometric", \
			Pvalue= self.hypergeometric(N, K, n, k)))
      else:
       	enrich.append(EnrichmentGFP(gfp_id=self._gfpid, Profile=profile, \
			ProfileID=int(self._profileid[i]), Test="hypergeometric", Pvalue= 1.0))
      i += 1

    if EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile, Test="hypergeometric").exists():
      enrichment = EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile).order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)
    else:
      EnrichmentGFP.objects.bulk_create(enrich, batch_size=500)

      enrichment = EnrichmentGFP.objects.filter(gfp_id=self._gfpid, Profile=profile, Test="hypergeometric").order_by("Pvalue")
      pvalues = []
      for enrich in enrichment:
        pvalues.append(enrich.Pvalue)

      pvaluesbh = self.BH(pvalues)
      i = 0
      for enrich in enrichment:
        enrich.BHPvalue = pvaluesbh[i]
        i += 1
      bulk_update(enrichment, batch_size=500)

  def hypergeometricown(self, N, K, n, k):
    """
      N= total number of genes in population
      K= number of GOA
      n= select a sample (top 50, bottom half, etc.)
      k= number of successes in the sample
    """
    return comb(K, k) * comb(N-K, n-k) / comb(N, n)
  def hypergeometric(self, N, K, n, k):
    """
      N= total number of genes in population
      K= number of GOA
      n= select a sample (top 50, bottom half, etc.)
      k= number of successes in the sample
    """
    return hypergeom(N, K, n).pmf(k)
  def hypergeometric_cdf(self, N, K, n, k):
    """
      N= total number of genes in population
      K= number of GOA
      n= select a sample (top 50, bottom half, etc.)
      k= number of successes in the sample
    """
    return 1 - hypergeom(N, K, n).cdf(k)

  def BH(self, pvalues):
    return multipletests(pvalues, method="fdr_bh")[1]

  def getOverlap(self, ProfileID, novel=False):
    #print "ProfileID", ProfileID
    #print self._profileid.shape
    #print self._profile.shape

    index = where(self._profileid == ProfileID)[0]
    #print "index: ", index

    function = asarray(self._profile[:, index]).ravel()
    #print "function: ", function.shape
    #function = function.ravel()
    #(gosize,) = where(self._profile[:, index] == 1)[0].shape
    (gosize,) = where(function == 1)[0].shape

    #print "novel: ", self._novel.shape

    overlap = zeros(function.shape)
    if novel:
      overlap = function + self._novel
    else:
      overlap = function + self._annotation

    indeces = where(overlap == 2)[0]
    (qoverlap,) = indeces.shape

    proteinover = self._proteinid[indeces]

    genes = []
    for protein in proteinover:
      enriched = NovelCandidates.objects.filter(ProteinID=protein, gfp=self._gfpid)
      genes.append(enriched[0].Symbol)
    return (gosize, qoverlap, genes)

  def MTCown(self, pvalues, correction_type = "Benjamini-Hochberg"):                
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    #from numpy import array, empty                                                                        
    pvalues = array(pvalues) 
    n = float(pvalues.shape[0])                                                                           
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":                                                                   
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":                                                            
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        for rank, vals in enumerate(values):                                                              
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue                                                            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * pvalue)                                                          
        for i in xrange(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                           
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues

class TestEnrichmentAnalysis(TestCase):
  def setUp(self):

    import django
    django.setup()
    #self.enrich = EnrichmentAnalysis(1, 131)
    self.enrich = EnrichmentAnalysis(18, 170)
    #self.enrich = EnrichmentAnalysis(9, 135)
  def tearDown(self):
    del self.enrich
  def testAEnrich(self):

    #self.enrich.loadProfile("GOID_Human_1_24_2015.mat", \
    #				"Gene2GO_Human_1_24_2015.mat")

    #self.enrich.run()
    pass

  def testBscRNASeqEnrich(self):
    #enrich = EnrichmentAnalysis(170, 18)
    #self.enrich.loadProfile("GOID_Human_1_24_2015.mat", \
    #				"Gene2GO_Human_1_24_2015.mat")
    #self.enrich.runscRNASeq(profile="GO")
    #self.enrich.runGFPscRNASeq(profile="GO")
    #del enrich

    #self.enrich.loadProfile("KEGGID_Human_1_24_2015.mat", \
    #				"KEGG_Human_1_24_2015.mat")
    #self.enrich.runGFPscRNASeq(profile="KEGG")
    pass


def suite():
        suite = makeSuite(TestEnrichmentAnalysis, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])

