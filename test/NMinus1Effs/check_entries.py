import sys
import os
import shutil
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram

ROOT.TH1.AddDirectory(False)

# dirs = [	
# 	'/afs/cern.ch/work/f/ferrico/private/Codice_ZPrime_8_NoTrigger/CMSSW_8_0_3_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/NMinus1Effs/data'
# 	'/afs/cern.ch/work/f/ferrico/private/Codice_ZPrime_8_NoTrigger/CMSSW_8_0_3_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/NMinus1Effs/mc'
# 		]

#file_to_check = '/afs/cern.ch/work/f/ferrico/private/Codice_ZPrime_8_NoTrigger/CMSSW_8_0_3_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/NMinus1Effs/mc/ana_nminus1_%s.root'
file_to_check = '/afs/cern.ch/work/f/ferrico/private/Codice_ZPrime_8_NoTrigger/CMSSW_8_0_3_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/mc/ana_datamc_%s.root'
histograms = ['DimuonMassVertexConstrained', 'DimuonMassVertexConstrained_bb', 'DimuonMassVertexConstrained_ne', 'DimuonMassVertexConstrained_pe' 
				]

#samples = ['dy50to120', 'dy120to200','dy200to400','dy400to800','dy800to1400','dy1400to2300','dy2300to3500','dy3500to4500','dy4500to6000','WWinclusive', 'WW200to600', 'WW600to1200', 'WW1200to2500', 'WW2500', 'WZ','ZZ', 'Wantitop','tW', 'Wjets', 'qcd80to120', 'qcd120to170', 'qcd170to300', 'qcd300to470', 'qcd470to600', 'qcd600to800','qcd800to1000','qcd1000to1400','qcd1400to1800','qcd1800to2400','qcd2400to3200','qcd3200']#,ttbar_pow, 'dyInclusive50
samples = ['dy50to120', 'dy120to200','dy200to400','dy400to800','dy800to1400','dy1400to2300','dy2300to3500','dy3500to4500','dy4500to6000', 'WW200to600', 'WW600to1200', 'WW1200to2500', 'WW2500','WZ','ZZ', 'Wantitop','tW', 'ttbar_pow']#,ttbar_pow, 'dyInclusive50, 'Wjets'

cuts = [
			'NoPt',
		    'NoDB',
		    'NoIso',
		    'NoTkLayers',
		    'NoPxHits',
		    'NoMuHits',
		    'NoMuMatch',
		    'NoVtxProb',
		    'NoB2B',
		    'NoDptPt',
		    'NoCosm',
		    'NoTrgMtch',
		    'NoNo',
		    'TiDB',
		    'TiIso',
		    'TiMuHits',
		    'TiMuMatch',
		    'TiPt',
		    'TiPxHits',
		    'TiTkLayers'
		]
		



#file = ROOT.TFile(file_to_check)	
# somma = 0
# for histo in histograms:
# 	totale = file.NoDB.Get(histo).Clone()
# 	print '%-35s --- %-7d' % (histo, totale.Integral())
# 	if histo == 'DimuonMassVertexConstrained':
# 		somma = somma - totale.Integral()
# 	else:
# 		somma = somma + totale.Integral()
# print '%63d' % somma
	
		
# i = 1 
# for sample in samples:
# # 	print sample
# 	print i
# 	somma = 0
# 	file = ROOT.TFile(file_to_check % sample)
# 	for histo in histograms:
# 		totale = file.NoVtxProb.Get(histo).Clone()
# 		print '%-15s --- %-35s --- %-7d' % (sample, histo, totale.Integral())
# 		if histo == 'DimuonMassVertexConstrained':
# 			somma = somma - totale.Integral()
# 		else:
# 			somma = somma + totale.Integral()
# 	print '%83d' % somma
# 	if somma != 0:
# 		print ' *************     ALLERT    ************'
# 	i = i+1
	
i = 1 
for sample in samples:
# 	print sample
	print i
	somma = 0
	file = ROOT.TFile(file_to_check % sample)
	totale = file.SimpleMuonsPlusMuonsMinusHistos.Get('DimuonMassVertexConstrained').Clone()
	bb = file.SimpleMuonsPlusMuonsMinusBarrelHistos.Get('DimuonMassVertexConstrained').Clone()
	ne = file.SimpleMuonsPlusMuonsMinusNegativeHistos.Get('DimuonMassVertexConstrained').Clone()
	pe = file.SimpleMuonsPlusMuonsMinusPositiveHistos.Get('DimuonMassVertexConstrained').Clone()
	print sample
	print 'Tot = %15d' % totale.Integral()
	print ' bb = %15d' % bb.Integral()
	print ' ne = %15d' % ne.Integral()

	print ' pe = %15d' % pe.Integral()
	somma = totale.Integral() - (bb.Integral() + ne.Integral() + pe.Integral())
	print '%25d' % somma
	if somma != 0:
		print '                                   *************    ALLERT: %15s    ************' % sample
	i = i+1


