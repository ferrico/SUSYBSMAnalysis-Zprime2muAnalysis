#!/usr/bin/env python

import os
import array
import numpy
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)

# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(111)

# variable = 'DimuonMassVertexConstrained'
# variable = 'DimuonMassVertexConstrained_bb'
# category = 'Barrel - Barrel'
# variable = 'DileptonMass_bb'
variable = 'LeptonPt'
# category = 'BE + EE'
# variable = 'DimuonMassVertexConstrained_be'

low = fitlow = 150
high = fithigh = 2500
high_for_data = 2500

bias = ['0p01', '0p05', '0p1', '0p15', '0p2', '0p5']

latex = ROOT.TLatex()

int_lumi = 36238.734
rescale_factor = 0.9329
rebin = 40
bin = array.array('d', [150, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 2500]);
use_non_dy = True

masses  = ['dy50to120', 'dy120to200', 'dy200to400', 'dy400to800', 'dy800to1400', 'dy1400to2300', 'dy2300to3500', 'dy3500to4500', 'dy4500to6000']
nevents = [2977600, 100000, 100000,  98400, 100000, 100000, 100000, 100000, 100000]
sigmas  = [  1975, 19.32, 2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7]
weights = [rescale_factor * int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]

hists_NoScale = []
hists_dir = '../DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/'
for m,w in zip(masses, weights):
	fn = 'ana_datamc_%s.root' % m
	fn = hists_dir + fn
	f = ROOT.TFile(fn)
	d = f.Our2016MuonsPlusMuonsMinusHistos
	h = d.Get(variable).Clone('%s' % m)
	h.Rebin(rebin)
	h.GetXaxis().SetRangeUser(low, high)
# 	print m,w
	h.Scale(w)
	h.Draw()
# 	ps.save('rawmass%s' % m)
	hists_NoScale.append(h)
 
if use_non_dy:
	from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
	non_dy_samples = [		WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, 
    						WZ,#WZ_skim,
							ZZ,#ZZ_skim,
    						WZ_ext, 
    						ZZ_ext,#ZZ_ext_skim,
    						Wantitop, tW, 
    						Wjets, 
    						ttbar_lep50to500, 
							ttbar_lep_500to800, 
							ttbar_lep_800to1200, 
							ttbar_lep_1200to1800, 
							ttbar_lep1800toInf,
							#     						qcd50to80, 
    						qcd80to120, qcd120to170, 
    						qcd170to300, 
    						qcd300to470, qcd470to600, qcd600to800, qcd800to1000, qcd1000to1400, 
    						qcd1400to1800, qcd1800to2400, qcd2400to3200, qcd3200, 
    						dyInclusive50
    						]
for sample in non_dy_samples:
       	fn = 'ana_datamc_%s.root' % sample.name
       	fn = hists_dir + fn
       	f = ROOT.TFile(fn)
       	d = f.Our2016MuonsPlusMuonsMinusHistos
       	w = sample.partial_weight * int_lumi
       	h = d.Get(variable).Clone('%s' % sample.name)
       	h.Rebin(rebin)
       	h.GetXaxis().SetRangeUser(low, high)
       	h.Scale(w)
       	h.Draw()
#        	ps.save('rawmass_%s' % sample.name)
#        	print sample.name, w
       	hists_NoScale.append(h)
  
htot_NoScale = hists_NoScale[0].Clone()
	
for j in xrange(1, len(hists_NoScale)):
	htot_NoScale.Add(hists_NoScale[j])
	
htot_NoScale.SetLineColor(ROOT.kBlack)
htot_NoScale.Draw()

for bia in bias:
	ps = plot_saver('plots/SCALE_' + bia + '/' + variable)
	print bia
	hists = []
	hists_dir_DY = '../DataMCSpectraComparison/mc/DY_kFactor_Scale_' + bia + '/'
# 	hists_dir = '../DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/'
	for m,w in zip(masses, weights):
		fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
# 		if bia == 'NoScale':
# 			fn = hists_dir + fn
# 		else:
		fn = hists_dir_DY + fn
		f = ROOT.TFile(fn)
		d = f.Our2016MuonsPlusMuonsMinusHistos
		h = d.Get(variable).Clone('%s' % m)
		h.Rebin(rebin)
		h.GetXaxis().SetRangeUser(low, high)
# 		print m,w
		h.Scale(w)
		h.Draw()
		ps.save('rawmass%s' % m)
		hists.append(h)
 
	if use_non_dy:
		from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
		non_dy_samples = [		WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, 
    						WZ,#WZ_skim,
							ZZ,#ZZ_skim,
    						WZ_ext, 
    						ZZ_ext,#ZZ_ext_skim,
    						Wantitop, tW, 
    						Wjets, 
    						ttbar_lep50to500, 
							ttbar_lep_500to800, 
							ttbar_lep_800to1200, 
							ttbar_lep_1200to1800, 
							ttbar_lep1800toInf,
							#     						qcd50to80, 
    						qcd80to120, qcd120to170, 
    						qcd170to300, 
    						qcd300to470, qcd470to600, qcd600to800, qcd800to1000, qcd1000to1400, 
    						qcd1400to1800, qcd1800to2400, qcd2400to3200, qcd3200, 
    						dyInclusive50
    						]
    	for sample in non_dy_samples:
        	fn = 'ana_datamc_%s.root' % sample.name
	        fn = hists_dir + fn
	        f = ROOT.TFile(fn)
	        d = f.Our2016MuonsPlusMuonsMinusHistos
	        w = sample.partial_weight * int_lumi
	        h = d.Get(variable).Clone('%s' % sample.name)
	        h.Rebin(rebin)
	        h.GetXaxis().SetRangeUser(low, high)
	        h.Scale(w)
	        h.Draw()
	        ps.save('rawmass_%s' % sample.name)
# 	        print sample.name, w
	        hists.append(h)
  
  
	htot = hists[0].Clone()
	
	for j in xrange(1, len(hists)):
		htot.Add(hists[j])

	print "ciola"
	
	htot.SetTitle('')
	htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
	htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
	
	htot_NoScale.SetLineColor(ROOT.kBlack)
	htot_NoScale.Draw()
	htot.Draw("same")
	htot.SetLineColor(ROOT.kRed)
	legend = ROOT.TLegend(0.6, 0.75, 0.9, 0.95)
	legend.AddEntry(htot_NoScale, "MC: DY NO bias", "l")
	legend.AddEntry(htot, "MC: DY with %s bias" % bia, "l")
	legend.Draw()
    
	def fit_it(lo, hi):
				
		print " ----------------------------------------------------------------------------------====================== INIZIO "
		
		htot_NoScale.Draw()
		htot.Draw("same")

		ps.c.Update()
		
		fcn = ROOT.TF1("fcn", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
				(x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
		fcn.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
		fcn.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
		fcn.SetLineColor(ROOT.kBlue)
		if 'LeptonPt' != variable:
			htot.Fit(fcn, 'SREMV')    
		s = htot.GetListOfFunctions().FindObject("stats")
		s.SetName("Const")
		s.SetX1NDC(0.73)
		s.SetY1NDC(0.69)
		s.SetY2NDC(0.94)
		s.SetOptStat(10)
		s.SetOptFit(11111)
		s.SetTextColor(ROOT.kBlue)
		
		htot.SetMinimum(10e-6)
		htot.SetMaximum(10e5)
		
		ps.c.Update()   
		ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
		
		htot.GetXaxis().SetRangeUser(600,1200)
		
		ps.save('mass%i_%i_ZOOM' % (600, 1200), pdf=True, pdf_log=True)
		ps.c.Update()
		
		xax = htot.GetXaxis()
		if 'LeptonPt' != variable:
			hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-MC)/MC', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
		else:
			hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(NoScale - Scale)/ NoScale', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
		
		for h in [hres]:
			h.SetMarkerStyle(2)
			h.SetStats(0)
			
		for i in xrange(1, hres.GetNbinsX()+1):
			xlo = xax.GetBinLowEdge(i)
			xhi = xax.GetBinLowEdge(i+1)
			if xlo >= lo and xhi <= hi:#and xhi < 500:
				if 'LeptonPt' == variable:
					res = htot_NoScale.GetBinContent(i) - htot.GetBinContent(i)
				else:
					res = fcn.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
	
				if htot.GetBinContent(i) > 0:
					if 'LeptonPt' == variable:
						hres.SetBinContent(i, res/htot_NoScale.GetBinContent(i))
						hres.SetBinError(i, htot_NoScale.GetBinError(i)/htot_NoScale.GetBinContent(i))
					else:
						hres.SetBinContent(i, res/htot.GetBinContent(i))
						hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))
		hres.SetMinimum(-1)
		hres.SetMaximum(1)
		hres.GetXaxis().SetRangeUser(lo, hi)
		hres.Draw('e')
		l1 = ROOT.TLine(lo, 0., hi,  0.)
		l1.Draw()
		
		ps.save('res_%i_%i' % (lo, hi), pdf=True, log=False)
		ps.c.Update()
		
		f_data = ROOT.TFile('/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data_OK/NoScale_YesEtaCut_Run2016MuonsOnly/ana_datamc_data.root')
		c_data = f_data.Our2016MuonsPlusMuonsMinusHistos
		data = c_data.Get(variable)
		data.Rebin(rebin)
		data.SetStats(0)
		data.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
		data.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
		data.SetLineColor(ROOT.kRed)
		htot.SetLineColor(ROOT.kBlue)
		htot.GetXaxis().SetRangeUser(lo, hi)
		htot.Draw()
		if 'LeptonPt' != variable:
			htot.Fit(fcn, 'SREMV')
		data.Draw("same")
		legend = ROOT.TLegend(0.3, 0.75, 0.6, 0.85)
		legend.AddEntry(data, "Data", "l")
		legend.AddEntry(htot, "MC", "l")
		legend.AddEntry(fcn, "Fit MC", "l")
		legend.Draw()
		htot.GetXaxis().SetRangeUser(lo, high_for_data)
		data.GetXaxis().SetRangeUser(lo, high_for_data)
		data.SetMaximum(10e5)
		data.SetMinimum(10e-6)
		
		ps.save('data_%i_%i' % (lo, hi), pdf_log=True)
		
		data.SetMinimum(10)
		data.GetXaxis().SetRangeUser(600,1200)
		ps.save('data%i_%i_ZOOM' % (600, 1200), pdf_log=True)
		
		data.GetXaxis().SetRangeUser(lo, high_for_data)
		if 'LeptonPt' == variable:
			hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(MC-data)/data', data.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(data.GetNbinsX()+1))
		else:
			hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-data)/data', data.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(data.GetNbinsX()+1))
		xax = data.GetXaxis()
		for h in [hres]:
			h.SetMarkerStyle(2)
			h.SetStats(0)
#        h.GetYaxis().SetLabelSize(0.02)
		for i in xrange(1, data.GetNbinsX()+1):
			xlo = xax.GetBinLowEdge(i)
			xhi = xax.GetBinLowEdge(i+1)
			if xlo >= lo and xhi <= hi:
				if 'LeptonPt' != variable:
					res = fcn.Integral(xlo, xhi)/(xhi-xlo) - data.GetBinContent(i)
				else:
					res = htot.GetBinContent(i) - data.GetBinContent(i)
# 				print " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %f" %res
				if data.GetBinContent(i) > 0:
					hres.SetBinContent(i, res/data.GetBinContent(i))
					hres.SetBinError(i, data.GetBinError(i)/data.GetBinContent(i))
				if data.GetBinContent(i) == 0:
					hres.SetBinContent(i, 1)
					hres.SetBinError(i, 0)
		hres.SetMinimum(-1)
		hres.SetMaximum( 1)
		hres.GetXaxis().SetRangeUser(lo, high_for_data)
		hres = hres.Rebin(12, "ciola", bin);
		ps.c.Update() 
		hres.Draw('e')
		l1 = ROOT.TLine(lo, 0., high_for_data,  0.)
		l1.Draw()
		
		ps.save('res_%i_%i_data' % (lo, hi), pdf=True, log=False)
		ps.c.Update() 
		
    
#    print fcn.GetProb(), fcn.GetParameter(2), fcn.GetParError(2)
#     c = r.GetCovarianceMatrix()
#     r.Print("V")
#     print " ----------------------------------------------------------------------------------====================== FINE "
#     print c

	l = [150] 
	for lo in l:
		fit_it(lo, high)
		print " ------------------------------------------------------------------------- %s" %bia

print variable, high
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
# draw_overlay(fsets)
