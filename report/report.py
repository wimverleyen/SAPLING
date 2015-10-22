import logging

from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm, cm, inch
from reportlab.lib.pagesizes import letter, landscape, portrait
from reportlab.platypus import Paragraph, Table, TableStyle, Image

from GFPs.models import GFP, Prediction, Performance, AggregatePerformance, AggregatePrediction, GeneSymbol, NovelCandidates, Enrichment, EnrichmentGFP, GO, Experiment, CoAnnotation
from enrich.enrichment import EnrichmentAnalysis

import os
import pickle
import csv
from scipy.io import loadmat, savemat
from scipy.stats import rankdata
import numpy as np

from sapling.settings import BASE_DIR

log = logging.getLogger(__name__)

FIG_DIR = os.path.join(BASE_DIR, "figures/")
REPORT_DIR = os.path.join(BASE_DIR, "report/")


class Report:
  def __init__(self, experimentID):

    self.width, self.height = letter
    self.experimentid = experimentID

    ProteinID="entrezGeneID_Human_1_24_2015.mat"
    path = BASE_DIR + "/data/" + ProteinID
    proteinid = loadmat(path, squeeze_me=False, mat_dtype=True)
    self._proteinid = proteinid['proteinid']
    self.gfpid = GFP.objects.filter(experiment_id=experimentID)

    self.map = {}

    mapper="Mapper_Symbol_Entrez.pkl"
    path = BASE_DIR + "/data/Mappers/" + mapper
    handle = open(path, 'r')
    symbolsentrez = pickle.load(handle)
    handle.close()

    self.style = ParagraphStyle('Body', fontName="Helvetica", fontSize=9, spaceBefore=24)

    self.tablestyle = TableStyle([ 
      ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#EDEDED')),
      ('FONTSIZE', (0, 0), (-1, -1), 7), ])

    i = 0
    for m in symbolsentrez.keys():
      self.map[int(symbolsentrez[m])] = m

    self.types = {}
    self.types["Co"] = "co-expression"
    self.types["PPI"] = "protein-protein interaction"
    self.types["SM"] = "semantic similarity"

    self.algorithms = {}
    self.algorithms["NV"] = "neighbor voting"
    self.algorithms["LR"] = "logistic regression"
    self.algorithms["RWR"] = "random walk with restart"
    self.algorithms["SGD"] = "support vector machines"
    self.algorithms["PA"] = "passive aggressive"
    self.algorithms["RP"] = "rankprop"

  def coord(self, x, y, unit=1):
    """
      # http://stackoverflow.com/questions/4726011/wrap-text-in-a-table-reportlab
      Helper class to help position flowables in Canvas objects
    """
    x, y = x * unit, self.height -  y * unit
    return x, y

  def addPageNumber(self, canvas):
    txt = "%s" % canvas.getPageNumber()
    canvas.setFont('Helvetica', 9, leading=None)
    canvas.drawRightString(575, 20, txt)

  def enrichCSV(self, profile="GO"):
    filename = REPORT_DIR + "Enrichment_GeneList_GO_" + str(self.experimentid) + ".csv"
    enrich = Enrichment.objects.filter(experiment_id=self.experimentid, Profile="GO").order_by("BHPvalue")[0:50]

    enriched = EnrichmentAnalysis(self.experimentid, self.gfpid[0])
    enriched.loadProfile("GOID_Human_1_24_2015.mat", \
			"Gene2GO_Human_1_24_2015.mat")

    lines = []
    header = "Identifier, Function, Enriched genes, GO term size, # overlap\n"
    lines.append(header)
    for function in enrich:
      line = ""
      line += str(function.ProfileID) + ", "

      if profile == "GO":
        go = GO.objects.filter(GOID=function.ProfileID)
        goname = go[0].Name
        line += goname + ", "
      
      (gosize, qoverlap, genes) = enriched.getOverlap(function.ProfileID)

      geneann = ""
      for gene in genes:
        geneann += gene + "; "
      line += geneann + ", "
      line += str(gosize) + ", "
      line += str(qoverlap) + "\n"

      lines.append(line)
      del line

    handler = open(filename, 'w')
    handler.writelines(lines)
    handler.close()

  def gfpCSV(self, gfpID):
    pass

  def analyzePerformance(self, AUROC, network, meanAUROC=0.0):
    analysis = ""
    if AUROC < 0.6:
      analysis += "The performance of this method is poor (AUROC &lt 0.6). "
    elif AUROC > 0.6 and AUROC < 0.7:
      analysis += "The performance of this method is moderate (0.6 &lt AUROC &lt 0.7). "
    elif AUROC > 0.7 and AUROC < 0.8:
      analysis += "The performance of this method is good (0.7 &lt AUROC &lt 0.8). "
    elif AUROC > 0.8 and AUROC < 0.9:
      analysis += "The performance of this method is very good (0.8 &lt AUROC &lt 0.9). "
    elif AUROC > 0.9:
      analysis += "The performance of this method is very good (0.9 &lt AUROC &lt 1.0). "
    if AUROC < meanAUROC:
      analysis += "Predicting GO terms with 20 - 1000 annotations with the same %s network data and the baseline algorithm (i.e., neighbor voting based upon the guilt-by-association (GBA) principle) performs in average than this method on the original gene list." % network
    else:
      analysis += "Novel candidate genes are predicted with higher performance compared to the mean performance of the GO terms with 20 - 1000 annotations with the same %s network data and the baseline algorithm (i.e., neighbor voting based upon the guilt-by-association (GBA) principle)." % network
    return analysis

  def generateReport(self, filename):
    # Read in data
    voffset = 65

    # Generate report
    c = canvas.Canvas(filename, pagesize=portrait(letter))
 
    c.setFont('Helvetica', 30, leading=None)
    c.drawCentredString(315, 740, "SAPLING report")

    experiment = Experiment.objects.filter(id=self.experimentid)
    experiment = experiment[0]

    ptext = "This SAPLING report provides an overview of the analysis of your gene list regarding your %s experiment. First (1), your gene list will be conventionally analyzed for enrichment with gene ontology (GO) terms by using the Fisher's exact test. Second (2), SAPLING will predict novel associated candidate genes with user-defined gene function prediction methods. These methods use field-wide biological data resources represented in a network, i.e., nodes represent genes and edges represent the strength of the relationship between these genes, and a machine learning algorithm. In this section, we list the top 20 novel predicted candidate genes, provide an overview of the 3-fold cross-validation held-out procedure for the performance evaluation of the gene function prediction method, list top 20 novel candidate genes that contribute most to the method performance, provide the network visualization of the genes in the original gene list and the novel predicted candidate genes, and we analyze the network connectivity of the genes in the original gene list. At the end of this section we provide the enrichment analysis of the extended gene list, i.e., the original gene list extended with the union of the top 20 novel candidate genes and the top 20 novel candidate genes that contribute most to the method performance. The enrichment with gene ontology terms will be tested by using two statistical tests: Mann-Whitney (for the ranked list) and Fisher's exact test (for set of top ranked genes). Finally (3), the predictions from all these methods will be aggregated by averaging their output scores. The performance of the aggregated predictor will be analyzed and compared with its methods of origin." % experiment.Name
    p = Paragraph(ptext, self.style)
    p.wrapOn(c, self.width-70, self.height)
    p.drawOn(c, 30, 550)

    c.setFont('Helvetica', 12, leading=None)
    c.drawString(100, 520, "Summary of gene function prediction methods") 

    data = [['Algorithm', 'Data', 'Performance (AUROC)', 'Page']]

    p = 3
    i = 1
    for gfp in self.gfpid:

      sample = []
      sample.append(self.algorithms[gfp.Algorithm])
      sample.append(gfp.Network)
      performance = Performance.objects.filter(gfp=gfp)
      perf = performance[0].AUROC
      label = "%0.3lf" % perf
      sample.append(label)
      sample.append(str(p))
      data.append(sample)

      p += 5
      i += 1
      if i == 22:
        break

    step = (15 + (i * 15))

    t = Table(data, colWidths=[3.5 * cm, 3.5 * cm, 3.5 * cm, 2.0 * cm], style=self.tablestyle)
    t.wrapOn(c, self.width-70, self.height)
    t.drawOn(c, 150, (500 - step))

    if len(self.gfpid) > 1:
      performance = AggregatePerformance.objects.filter(experiment=self.experimentid)

      c.setFont('Helvetica', 12, leading=None)
      c.drawString(100, (500 - step - 25), "Aggregated gene function prediction method performance: %0.3lf" % performance[0].AUROC) 

    self.addPageNumber(c)
    c.showPage()

    c.setFont('Helvetica', 16, leading=None)
    c.drawString(50, 725, "1. GO Enrichment analysis")
    c.bookmarkPage("P1")
    c.addOutlineEntry("1. GO enrichment", "P1")

    c.setFont('Helvetica', 12, leading=None)
    c.drawCentredString(145, 650, "A. Gene ontology (hypergeometric)")
    c.bookmarkPage("P2")
    c.addOutlineEntry("A. Gene ontology (hypergeometric)", "P2", level=1)

    ptext = "During an enrichment analysis, the overlap between the original gene list and a 'known' gene list of a biological related function is quantified. Here, the enrichment analysis is performed on the gene onology (GO) terms with 20 - 1000 genes. All analysis is performed after propagation of GOA annotations through the ontology."
    p = Paragraph(ptext, self.style)
    p.wrapOn(c, self.width-70, self.height)
    p.drawOn(c, 30, 680)

    ptext = "For the 10 most enriched GO terms, the Benjamin-Hochberg corrected p-values of the hypergeometric statistical test are listed. This is followed by the genes contributing to enrichment, GO term size (i.e., the number of genes annotated by a given GO term), and the number of genes on the list with the corresponding GO terms (i.e., # overlap)."

    p = Paragraph(ptext, self.style)
    p.wrapOn(c, self.width-70, self.height)
    p.drawOn(c, 30, 600)

    data = [['Identifier', 'Function', 'BH corrected p-value'], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', '']]

    
    self.enrichCSV()

    enrich = Enrichment.objects.filter(experiment_id=self.experimentid, Profile="GO").order_by("BHPvalue")[0:10]

    j = 1
    for function in enrich:
      profile = function.Profile
      profileid = function.ProfileID
      data[j][0] = str(profileid)

      if profile == "GO":
        go = GO.objects.filter(GOID=profileid)
        goname = go[0].Name
        data[j][1] = goname[:80]

      pvalue = "%.2e" % function.BHPvalue
      data[j][2] = pvalue
      j += 1

    t = Table(data, colWidths=[1.5 * cm, 13.5 * cm, 3.0 * cm], style=self.tablestyle)
    t.wrapOn(c, self.width-70, self.height)
    t.drawOn(c, 50, 350)

    enriched = EnrichmentAnalysis(self.experimentid, self.gfpid[0])
    enriched.loadProfile("GOID_Human_1_24_2015.mat", \
			"Gene2GO_Human_1_24_2015.mat")

    data = [['Identifier', 'Enriched genes', 'GO term size', '# overlap'], \
            ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
            ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
            ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
            ['', '', '', '']]

    j = 1
    for function in enrich:
      data[j][0] = str(function.ProfileID)

      (gosize, qoverlap, genes) = enriched.getOverlap(function.ProfileID)
      genesymbols = ""
      k = 0
      for symbol in genes:
        if k < len(genes)-1:
          genesymbols += symbol + ", "
        else:
          genesymbols += symbol
        k += 1
      data[j][1] = genesymbols[:75]
      data[j][2] = str(gosize)
      data[j][3] = str(qoverlap)
      j += 1

    t = Table(data, colWidths=[1.5 * cm, 13 * cm, 2.0 * cm, 1.5 * cm], style=self.tablestyle)
    t.wrapOn(c, self.width-70, self.height)
    t.drawOn(c, 50, 100)

    self.addPageNumber(c)
    c.showPage()
  
    i = 0 
    gfpperf = []
    for gfp in self.gfpid:
      if i == 0:
        c.setFont('Helvetica', 16, leading=None)
        c.drawString(50, 725, "2. Gene function prediction")
        label = "P" + str((i*12) + 3 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("2. Gene function prediction", label)

        c.setFont('Helvetica', 12, leading=None)
        c.drawString(50, 700, "Method %d: %s network with %s algorithm" % (i+1, gfp.Network, \
					self.algorithms[gfp.Algorithm]))
        label = "P" + str((i*12) + 4 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("Method %d: %s network with %s algorithm" % (i+1, gfp.Network, \
					self.algorithms[gfp.Algorithm]), label, level=1)

        ptext = "In the following table, the top 20 novel associated genes predicted by the gene function prediction method composed of %s network data and %s algorithm are provided. For each predicted gene, the standardized rank (i.e., rank/max_rank) of the output score is reported." % (gfp.Network, \
			self.algorithms[gfp.Algorithm])

        p = Paragraph(ptext, self.style)
        p.wrapOn(c, self.width-70, self.height)
        p.drawOn(c, 30, 650)

      else:
        c.setFont('Helvetica', 12, leading=None)
        c.drawString(50, 725, "Method %d: %s network with %s algorithm" % (i+1, gfp.Network, \
					self.algorithms[gfp.Algorithm]))
        label = "P" + str((i*12) + 5 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("Method %d: %s network with %s algorithm" % (i+1, gfp.Network, \
					self.algorithms[gfp.Algorithm]), label, level=1)

        ptext = "In the following table, the top 20 novel associated genes predicted by the gene function prediction method composed of %s network data and %s algorithm are provided. For each predicted gene, the standardized rank (i.e., rank/max_rank) of the output score of the gene function prediction method is reported." % (gfp.Network, \
			self.algorithms[gfp.Algorithm])

        p = Paragraph(ptext, self.style)
        p.wrapOn(c, self.width-70, self.height)
        p.drawOn(c, 30, 660)

      data = [['Symbol', 'Score', 'Symbol', 'Score'], ['', '', '', ''], ['', '', '', ''], \
		['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
		['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
		['', '', '', ''], ['', '', '', '']]

      performance = Performance.objects.filter(gfp=gfp)
      predictions = NovelCandidates.objects.filter(gfp=gfp, Annotation=0, CV=-1).order_by('Score')[0:20]

      gfpperf.append(performance[0].AUROC)

      j = 1
      for prediction in predictions:

        symbol = prediction.Symbol
	rank = "%.5lf" % prediction.Rank

        if j < 11:
          data[j][0] = symbol
          data[j][1] = rank
        else:
          data[j-10][2] = symbol
          data[j-10][3] = rank

        j += 1

      t = Table(data, colWidths=[2.0 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm], style=self.tablestyle)
      t.wrapOn(c, self.width-70, self.height)
      t.drawOn(c, 150, 440)

      c.setFont('Helvetica', 12, leading=None)
      c.drawString(80, 410, "a. Performance evaluation")
      label = "P" + str((i*12) + 6 + 2)
      c.bookmarkPage(label)
      c.addOutlineEntry("a. Performance evaluation", label, level=2)

      if gfp.Type == 'Co':
	      path = BASE_DIR + '/data/' + gfp.Type + '/' + gfp.Network + '/GO_Performance_1_24_2015.mat'
      else:
	      path = BASE_DIR + '/data/' + gfp.Type + '/' + gfp.Network + '/GO_Performance_1_24_2015.mat'
      goperf = loadmat(path, squeeze_me=False, mat_dtype=True)
      goperf = goperf['perf']

      (x, ) = np.where(goperf > performance[0].AUROC)[0].shape

      p = x/float(109)

      roc = FIG_DIR + 'ROC_Curve_'+ str(gfp.id) +'.png'
      c.drawImage(roc, 50, -30, width=533, height=366)

      ptext = "The receiver operating characteristic (ROC) curve is used to evaluate the predictions of the gene function prediction method. Typically, the area under the ROC curve (AUROC) is applied as a metric to evaluate performance. The gene function prediction method is analyzed with a 3-fold cross-validation held-out procedure. For each fold, a randomly selected set of annotated genes are held-out and the method will be evaluated if the set of these held-out genes can be predicted to be associated to the reminder of the original gene list (mean AUROC over the 3 folds = %0.3lf). The left panel shows the mean ROC curve from the 3-fold cross-validation held-out procedure. The performance of the gene function prediction method on the original gene list is benchmarked by its performance for the gene ontology (GO) terms. The right panel shows the AUROC performance distribution with mean = %0.3lf and std = %0.3lf (empirical p ~ %0.3lf). " % (performance[0].AUROC, np.mean(goperf), np.std(goperf), p)
 
      analysis = self.analyzePerformance(performance[0].AUROC, gfp.Network, meanAUROC=np.mean(goperf))
      ptext += analysis

      p = Paragraph(ptext, self.style)
      p.wrapOn(c, self.width-70, self.height)
      p.drawOn(c, 30, 270)

      self.addPageNumber(c)
      c.showPage()

      c.setFont('Helvetica', 12, leading=None)
      c.drawString(80, 725, "b. Highest ranked genes during 3-fold cross-validation")
      label = "P" + str((i*12) + 7 + 2)
      c.bookmarkPage(label)
      c.addOutlineEntry("b. Highest ranked genes during 3-fold cross-validation", label, level=2)

      ptext = "During each fold of the 3-fold cross validation held-out procedure, we evaluate if the method is able to predict the held-out genes related to that fold based on the reminder of the original gene list. The predictions for each fold indicate which unannotated genes contribute most to the method performance. Therefore, we report an alternative list of top 20 novel candidate genes based on their rank during the cross-validation procedure."
      p = Paragraph(ptext, self.style)
      p.wrapOn(c, self.width-70, self.height)
      p.drawOn(c, 30, 660)

      data = [['Symbol', 'Score', 'Symbol', 'Score'], ['', '', '', ''], ['', '', '', ''], \
		['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
		['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
		['', '', '', ''], ['', '', '', '']]

      predictions = NovelCandidates.objects.filter(gfp=gfp, Annotation=0, CV=3).order_by('Score')[0:20]

      j = 1
      for prediction in predictions:

        symbol = prediction.Symbol
	rank = "%.5lf" % prediction.Rank

        if j < 11:
          data[j][0] = symbol
          data[j][1] = rank
        else:
          data[j-10][2] = symbol
          data[j-10][3] = rank

        j += 1

      t = Table(data, colWidths=[2.0 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm], style=self.tablestyle)
      t.wrapOn(c, self.width-70, self.height)
      t.drawOn(c, 150, 450)

      net = FIG_DIR + 'Network_'+ str(gfp.id) +'.png'
      c.drawImage(net, 25, 0, width=533, height=366)

      c.setFont('Helvetica', 12, leading=None)
      c.drawString(80, 425, "c. Network visualization")
      label = "P" + str((i*12) + 8 + 2)
      c.bookmarkPage(label)
      c.addOutlineEntry("c. Network visualization", label, level=2)

      ptext = "The gene function predictions can be shown within the network data from which they were derived. We visualize the interactions between the members of the original gene list (i.e, 'seed genes'; blue nodes) and the predicted top 20 novel candidate genes (red nodes). The direct interactions between the members of the original gene list and the novel candidate genes are colored in green. Further, we show the indirect interactions and the corresponding genes between the 'seed genes' and the novel candidates in gray."

      p = Paragraph(ptext, self.style)
      p.wrapOn(c, self.width-70, self.height)
      p.drawOn(c, 30, 350)

      self.addPageNumber(c)
      c.showPage()

      net = FIG_DIR + 'Node_Degree_'+ str(gfp.id) +'.png'
      c.drawImage(net, 25, 250, width=533, height=366)

      c.setFont('Helvetica', 12, leading=None)
      c.drawString(80, 725, "d. Connectivity")
      label = "P" + str((i*12) + 9 + 2)
      c.bookmarkPage(label)
      c.addOutlineEntry("d. Connectivity", label, level=2)

      c.setFont('Helvetica', 10, leading=None)
      ptext = "We use the node degree, i.e., the number of interactions of a gene in the network, to analyze the connectivity of the genes in the original gene list. The node degree and size of the gene set combine to provide us with a gene-specific null distribution of expected connectivities to other genes provided as training data. If the set of genes were preferentially linked, they will generally also have been predictable in cross-validation. We compare the expected connectivity of the genes in the original gene list in the complete %s network with their actual values in the subnetwork of the genes in the original gene list. In the left panel, we can indicate if the connectivity between genes in the original gene list is higher than expected from the complete network; the gray line draws the expected node degree, i.e., the ratio of the number of genes in the original gene list and the total number of genes in the network. The actual and predicted connectivity are compared. In the right panel, the probability density function for the residuals of the expected connectivity is shown. More positive values in the probability density function result in better method performance; higher interconnectivity between genes in the gene list makes it more likely that held-out genes will be predicted from the reminder of the gene list for each fold during the cross-validation procedure." % gfp.Network

      p = Paragraph(ptext, self.style)
      p.wrapOn(c, self.width-70, self.height)
      p.drawOn(c, 30, 570)

      self.addPageNumber(c)
      c.showPage()
      
      ## NEW
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

      if Type == 'sc RNA Seq':
        c.setFont('Helvetica', 12, leading=None)
        c.drawString(80, 725, "e. GO enrichment analysis")
        label = "P" + str((i*12) + 10 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("e. GO enrichment analysis", label, level=2)

        ptext = "After expanding our original gene list with predicted novel candidate genes, we can again perform a GO enrichment analysis on the extended gene list. As a reminder, we extend the original gene list with the union of the top 20 novel candidate genes predicted by the gene function prediction method and the top 20 novel candidate genes highest ranked during the 3-fold cross-validation held-out procedure."
        p = Paragraph(ptext, self.style)
        p.wrapOn(c, self.width-70, self.height)
        p.drawOn(c, 30, 650)


        c.setFont('Helvetica', 12, leading=None)
        c.drawString(100, 630, "1) Gene ontology (hypergeometric)")
        label = "P" + str((i*12) + 12 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("1) Gene ontology (hypergeometric)", label, level=3)

        ptext = "For 10 most enriched GO terms, the Benjamin-Hochberg corrected p-values of the hypergeometric statistical test are listed. This is followed by the enriched genes, GO term size (i.e., the number of genes annotated by a given GO term), and the number of enriched genes of the corresponding GO terms (i.e., # overlap) for the enriched GO terms. The total number of genes for this enrichment analysis is equal to the number of genes with a node degree bigger than zero."

        p = Paragraph(ptext, self.style)
        p.wrapOn(c, self.width-70, self.height)
        p.drawOn(c, 30, 560)

        data = [['Identifier', 'Function', 'BH corrected p-value'], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', '']]

        enrich = EnrichmentGFP.objects.filter(gfp_id=gfp, Test="hypergeometric", Profile="GO").order_by("BHPvalue")[0:10]

        j = 1
        for function in enrich:
          profile = function.Profile
          profileid = function.ProfileID
          data[j][0] = str(profileid)

          if profile == "GO":
            go = GO.objects.filter(GOID=profileid)
            data[j][1] = go[0].Name

          pvalue = "%.2e" % function.BHPvalue
          data[j][2] = pvalue
          j += 1

        t = Table(data, colWidths=[1.5 * cm, 13.5 * cm, 3.0 * cm], style=self.tablestyle)
        t.wrapOn(c, self.width-70, self.height)
        t.drawOn(c, 50, 320)

        data = [['Identifier', 'Enriched genes', 'GO term size', '# overlap'], \
              ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
              ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
              ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
              ['', '', '', '']]

        j = 1
        for function in enrich:
          data[j][0] = str(function.ProfileID)

          (gosize, qoverlap, genes) = enriched.getOverlap(function.ProfileID, novel=True)
          genesymbols = ""
          k = 0
          for symbol in genes:
            if k < len(genes)-1:
              genesymbols += symbol + ", "
            else:
              genesymbols += symbol
            k += 1
          data[j][1] = genesymbols[:75]
          data[j][2] = str(gosize)
          data[j][3] = str(qoverlap)
          j += 1

        t = Table(data, colWidths=[1.5 * cm, 13 * cm, 2.0 * cm, 1.5 * cm], style=self.tablestyle)
        t.wrapOn(c, self.width-70, self.height)
        t.drawOn(c, 50, 80)
 
        self.addPageNumber(c)
        c.showPage()

      else:

        c.setFont('Helvetica', 12, leading=None)
        c.drawString(80, 725, "e. GO enrichment analysis")
        label = "P" + str((i*12) + 10 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("e. GO enrichment analysis", label, level=2)

        ptext = "After expanding our original gene list with predicted novel candidate genes, we can again perform a GO enrichment analysis on the extended gene list. As a reminder, we extend the original gene list with the union of the top 20 novel candidate genes predicted by the gene function prediction method and the top 20 novel candidate genes highest ranked during the 3-fold cross-validation held-out procedure. Since a gene function prediction method provides a continuous score for each gene, we can also use the Mann-Whitney statistical test for the enrichment analysis."
        p = Paragraph(ptext, self.style)
        p.wrapOn(c, self.width-70, self.height)
        p.drawOn(c, 30, 650)

        c.setFont('Helvetica', 12, leading=None)
        c.drawString(100, 630, "1) Gene ontology (Mann-Whitney)")
        label = "P" + str((i*12) + 11 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("1) Gene ontology (Mann-Whitney)", label, level=3)

        ptext = "The significance of the enrichment of the output scores of the gene function prediction method for each GO term is defined by the corresponding Benjamin-Hochberg corrected p-values from the Mann-Whitney test. Here, we present the 10 most enriched GO terms."
        p = Paragraph(ptext, self.style)
        p.wrapOn(c, self.width-70, self.height)
        p.drawOn(c, 30, 590)

        data = [['Identifier', 'Function', 'BH corrected p-value'], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', '']] 

        enrich = EnrichmentGFP.objects.filter(gfp_id=gfp, Test="Mann-Whitney", Profile="GO").order_by("BHPvalue")[0:10]

        j = 1
        for function in enrich:
          profile = function.Profile
          profileid = function.ProfileID
          data[j][0] = str(profileid)

          if profile == "GO":
            go = GO.objects.filter(GOID=profileid)
            goname = go[0].Name
            data[j][1] = goname[:80]

          pvalue = "%.2e" % function.BHPvalue
          data[j][2] = pvalue
          j += 1

        t = Table(data, colWidths=[1.5 * cm, 13.5 * cm, 3.0 * cm], style=self.tablestyle)
        t.wrapOn(c, self.width-70, self.height)
        t.drawOn(c, 50, 360)

        enriched = EnrichmentAnalysis(self.experimentid, gfp)
        enriched.loadProfile("GOID_Human_1_24_2015.mat", \
     			"Gene2GO_Human_1_24_2015.mat")

        self.addPageNumber(c)
        c.showPage()

        c.setFont('Helvetica', 12, leading=None)
        c.drawString(100, 725, "2) Gene ontology (hypergeometric)")
        label = "P" + str((i*12) + 12 + 2)
        c.bookmarkPage(label)
        c.addOutlineEntry("2) Gene ontology (hypergeometric)", label, level=3)

        ptext = "For 10 most enriched GO terms, the Benjamin-Hochberg corrected p-values of the hypergeometric statistical test are listed. This is followed by the enriched genes, GO term size (i.e., the number of genes annotated by a given GO term), and the number of enriched genes of the corresponding GO terms (i.e., # overlap) for the enriched GO terms."

        p = Paragraph(ptext, self.style)
        p.wrapOn(c, self.width-70, self.height)
        p.drawOn(c, 30, 680)

        data = [['Identifier', 'Function', 'BH corrected p-value'], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', ''], ['', '', ''], ['', '', ''], \
            ['', '', '']]

        enrich = EnrichmentGFP.objects.filter(gfp_id=gfp, Test="hypergeometric", Profile="GO").order_by("BHPvalue")[0:10]

        j = 1
        for function in enrich:
          profile = function.Profile
          profileid = function.ProfileID
          data[j][0] = str(profileid)

          if profile == "GO":
            go = GO.objects.filter(GOID=profileid)
            data[j][1] = go[0].Name

          pvalue = "%.2e" % function.BHPvalue
          data[j][2] = pvalue
          j += 1

        t = Table(data, colWidths=[1.5 * cm, 13.5 * cm, 3.0 * cm], style=self.tablestyle)
        t.wrapOn(c, self.width-70, self.height)
        t.drawOn(c, 50, 400)

        data = [['Identifier', 'Enriched genes', 'GO term size', '# overlap'], \
              ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
              ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
              ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
              ['', '', '', '']]

        j = 1
        for function in enrich:
          data[j][0] = str(function.ProfileID)

          (gosize, qoverlap, genes) = enriched.getOverlap(function.ProfileID, novel=True)
          genesymbols = ""
          k = 0
          for symbol in genes:
            if k < len(genes)-1:
              genesymbols += symbol + ", "
            else:
              genesymbols += symbol
            k += 1
          data[j][1] = genesymbols[:75]
          data[j][2] = str(gosize)
          data[j][3] = str(qoverlap)
          j += 1

        t = Table(data, colWidths=[1.5 * cm, 13 * cm, 2.0 * cm, 1.5 * cm], style=self.tablestyle)
        t.wrapOn(c, self.width-70, self.height)
        t.drawOn(c, 50, 150)
 
        self.addPageNumber(c)
        c.showPage()

      i += 1

    if len(self.gfpid) > 1:

      performance = AggregatePerformance.objects.filter(experiment=self.experimentid)
      predictions = AggregatePrediction.objects.filter(experiment=self.experimentid, Annotation=0, CV=-1).order_by('Score')[0:20]

      c.setFont('Helvetica', 16, leading=None)
      c.drawString(50, 725, "3. Aggregated gene function prediction")
      label = "P" + str((i*12) + 13 + 2)
      c.bookmarkPage(label)
      c.addOutlineEntry("3. Aggregated gene function prediction", label)

      analysis = ""
      analysis2 = ""
      if len(gfpperf) == 1:
        analysis += "The is no aggregation applied since there is only one gene function prediction configured in your experiment. "
      else:
        increase = 0.0
        if np.max(gfpperf) > 0.5:
          increase = (performance[0].AUROC - 0.5)/float((np.max(gfpperf) - 0.5))
        xx = 0.0
        if increase < 1.0:
          xx = (1.0 - increase) * 100
          analysis += "The aggregate predictor has an AUROC decrease below the maximum of the methods of %0.3lf %%. " % xx
        else:
          xx = (increase - 1.0) * 100
          analysis += "The aggregate predictor has an AUROC increase above the maximum of the methods of %0.3lf %%. " % xx
      
        increase = 0.0
        if np.max(gfpperf) > 0.5:
          increase = (performance[0].AUROC - 0.5)/float((np.mean(gfpperf) - 0.5))
        xx = 0.0
        if increase < 1.0:
          xx = (1.0 - increase) * 100
          analysis2 += "and an AUROC decrease below the mean of the methods of %0.3lf %%." % xx
        else:
          xx = (increase - 1.0) * 100
          analysis2 += "and an AUROC increase above the mean of the methods of %0.3lf %%." % xx

      ptext = "Finally, the 20 highest ranked genes derived by averaging the output scores from all gene function prediction methods are listed. Aggregation of the output scores results in an AUROC of %0.3lf. " % (performance[0].AUROC)

      ptext += analysis
      ptext += analysis2

      p = Paragraph(ptext, self.style)
      p.wrapOn(c, self.width-70, self.height)
      p.drawOn(c, 30, 665)

      data = [['Symbol', 'Score', 'Symbol', 'Score'], ['', '', '', ''], ['', '', '', ''], \
            ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
            ['', '', '', ''], ['', '', '', ''], ['', '', '', ''], \
            ['', '', '', ''], ['', '', '', '']]

      j = 1
      for prediction in predictions:
        symbols = GeneSymbol.objects.filter(ProteinID=str(prediction.ProteinID))
        symbol = ""

        for symbol in symbols:
          symbol = symbol.Symbol
          break

        rank = "%.5lf" % prediction.Rank

        if j < 11:
          data[j][0] = symbol
          data[j][1] = rank
        else:
          data[j-10][2] = symbol
          data[j-10][3] = rank

        j += 1

      t = Table(data, colWidths=[2.0 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm], style=self.tablestyle)
      t.wrapOn(c, self.width-70, self.height)
      t.drawOn(c, 150, 420)

      roc = FIG_DIR + 'ROC_Curve_Aggregate_'+ str(self.experimentid) +'.png'
      c.drawImage(roc, 150, 100, width=300, height=207)

      self.addPageNumber(c)

    c.save()
