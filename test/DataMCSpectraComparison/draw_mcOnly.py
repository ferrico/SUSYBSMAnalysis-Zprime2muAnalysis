# Testing GIT
#!/usr/bin/env python
import sys
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import overall_prescale

from pprint import pprint
from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *

import sys, os, glob
from collections import defaultdict
from pprint import pprint
from optparse import OptionParser

ROOT.gStyle.SetOptStat(0)


mumu_dati_new = '/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data/Run2016MuonsOnly/ana_datamc_data.root'
mumu_dati_old = '/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data/NoScale_YesEtaCut_Run2016MuonsOnly/ana_datamc_data.root'

histogram = ['DimuonMassVertexConstrained','DimuonMassVertexConstrained_bb', 'DimuonMassVertexConstrained_be', 'LeptonPt', 'DileptonMass', 'DileptonMass_bb', 'DileptonMass_be']
rebin = 40


for mumu_histogram in histogram:
	f_dati_new = ROOT.TFile(mumu_dati_new)
	dati_new = f_dati_new.Our2016MuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
	dati_new.Rebin(rebin)
	dati_new.SetMarkerColor(1)
	dati_new.SetLineColor(1)
# 	dati_new.SetLineWidth(2)
	
	f_dati_old = ROOT.TFile(mumu_dati_old)
	dati_old = f_dati_old.Our2016MuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
	dati_old.Rebin(rebin)
	dati_old.SetMarkerColor(2)
	dati_old.SetLineColor(2)
# 	dati_old.SetLineWidth(2)
	
	canvas = ROOT.TCanvas(mumu_histogram, mumu_histogram)#, 210,45,900,600)
	pad1 = ROOT.TPad("pad1", "pad1", 0, 0.2, 1, 1.0)
	pad1.SetBottomMargin(0)
	pad1.SetGridx()
	pad1.Draw()
	pad1.cd() 
	pad1.SetLogy()
	
	dati_old.Draw()
	dati_new.Draw("same")
	
	
	Xmin = 0
	Xmax = 3000

	dati_old.SetMaximum(10e5)
	dati_new.SetMaximum(10e5)
	dati_old.SetMinimum(10e-2)
	dati_new.SetMinimum(10e-2)
	
 	dati_old.GetXaxis().SetRangeUser(Xmin, Xmax)
 	dati_new.GetXaxis().SetRangeUser(Xmin, Xmax)
 	
 	dati_old.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV}")
 	dati_old.GetYaxis().SetTitle("event / %s GeV" % rebin)
 	dati_new.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV}")
 	dati_new.GetYaxis().SetTitle("event / %s GeV" % rebin)
 	
#  	legend = ROOT.TLegend(0.6,0.85,0.9,0.9)
#  	legend.AddEntry(dati_old,"No scale", "l")	
#  	legend.AddEntry(dati_new, "Scale APPLIED: NEW", "l")
#  	legend.Draw()
 	
#  	old = dati_old.GetListOfFunctions().FindObject("stats")
#  	new = dati_new.GetListOfFunctions().FindObject("stats")
#  	
#  	old.SetX1NDC(0.73)
#  	old.SetY1NDC(0.69)
#  	old.SetY2NDC(0.94)
#  	new.SetX1NDC(0.73)
#  	new.SetY1NDC(0.44)
#  	new.SetY2NDC(0.69)
#  	new.Draw()
#  	old.Draw()
 	
	canvas.Update()
	canvas.Print('./Comparison_data/%s.png' % mumu_histogram)
	