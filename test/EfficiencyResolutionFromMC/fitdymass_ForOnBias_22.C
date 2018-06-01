#!/usr/bin/env python
from math import sqrt
import os
import array
import numpy
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)


# variable = 'DimuonMassVertexConstrained'
# variable = 'DimuonMassVertexConstrained_bb'
# category = 'Barrel - Barrel'
# variable = 'DileptonMass_bb'
variable = 'LeptonPt'
# category = 'BE + EE'
# variable = 'DimuonMassVertexConstrained_be'


if variable != 'LeptonPt':
	ROOT.gStyle.SetOptFit(111)
else:
	ROOT.gStyle.SetOptFit(0)
	ROOT.gStyle.SetOptStat(0)


low = fitlow = 150
high = fithigh = 2500	
high_for_data = 2500

if 'LeptonPt' == variable:
	high = high_for_data = fithigh = 1900

bias = ['0p01', '0p05', '0p1', '0p15', '0p2', '0p5']

# latex = ROOT.TLatex()

int_lumi = 36238.734

if '_bb' in variable:
	rescale_factor = 0.9880
elif '_be' in variable:
	rescale_factor = 0.9625
else:
	rescale_factor = 0.9714

if 'LeptonPt' == variable:
	rebin = 100
else:
	rebin = 40


bin = array.array('d', [150, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 2500]);
use_non_dy = True


masses  = ['dy50to120', 'dy120to200', 'dy200to400', 'dy400to800', 'dy800to1400', 'dy1400to2300', 'dy2300to3500', 'dy3500to4500', 'dy4500to6000']
nevents = [2977600, 100000, 100000,  98400, 100000, 100000, 100000, 100000, 100000]
sigmas  = [  1975, 19.32, 2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7]
weights = [rescale_factor * int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]


bias_label = ROOT.TPaveText(0.3, 0.85, 0.6, 0.95, "NDC");
bias_label.SetTextAlign(12)
bias_label.SetTextFont(42)
bias_label.SetTextSize(0.04)
bias_label.SetFillStyle(0)
bias_label.SetBorderSize(0)

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
       	w = sample.partial_weight * int_lumi * rescale_factor
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
	
# htot_NoScale.SetLineColor(ROOT.kBlack)
# htot_NoScale.Draw()

for bia in bias:

	label = bia.replace('p', '.')
	
	bias_label = ROOT.TPaveText(0.3, 0.85, 0.6, 0.95, "NDC");
	bias_label.SetTextAlign(12)
	bias_label.SetTextFont(42)
	bias_label.SetTextSize(0.04)
	bias_label.SetFillStyle(0)
	bias_label.SetBorderSize(0)
	bias_label.AddText("Scale bias = %s" %label)

# 	ps = plot_saver('plots/SCALE_' + bia + '/' + variable)

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
# 		ps.save('rawmass%s' % m)
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
# 	        ps.save('rawmass_%s' % sample.name)
# 	        print sample.name, w
	        hists.append(h)
  
  
	htot = hists[0].Clone()
	
	for j in xrange(1, len(hists)):
		htot.Add(hists[j])

	

# 	bias_label.Draw()
# 	legend.Draw()
    
	def fit_it(lo, hi):		
		
		
		canvas = ROOT.TCanvas(variable, variable, 300, 60, 1000, 650);
		pad1 = ROOT.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
		pad1.SetBottomMargin(0)
		pad1.SetGridx()
		pad1.SetLogy()
		pad1.Draw()
		pad1.cd() 
		
		htot_NoScale.SetLineColor(ROOT.kBlack)
		htot.SetLineColor(ROOT.kRed)
		legend = ROOT.TLegend(0.3, 0.75, 0.6, 0.85)
		legend.AddEntry(htot_NoScale, "MC: DY NO bias", "l")
		legend.AddEntry(htot, "MC: DY with %s scale bias" % label, "l")
		htot_NoScale.Draw()
		htot.Draw("same")
		htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
		if variable == 'LeptonPt':
			htot.GetXaxis().SetTitle("p_{T} [GeV]")
			htot_NoScale.GetXaxis().SetTitle("p_{T} [GeV]")
		htot.SetTitle('')
		htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) 
		htot_NoScale.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) 
		htot.GetXaxis().SetRangeUser(lo, hi)
		htot_NoScale.GetXaxis().SetRangeUser(lo, hi)
		legend.Draw()

		ps.c.Update()
		print " aaaaaa %s" %bia
		
		fcn = ROOT.TF1("fcn", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
				(x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
		fcn.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
		fcn.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
		fcn.SetLineColor(ROOT.kBlue)
# 		if 'LeptonPt' != variable:
# 			htot.Fit(fcn, 'SREMV')    
# 			s = htot.GetListOfFunctions().FindObject("stats")
# 			s.SetName("Const")
# 			s.SetX1NDC(0.73)
# 			s.SetY1NDC(0.69)
# 			s.SetY2NDC(0.94)
# 			s.SetOptStat(10)
# 			s.SetOptFit(11111)
# 			s.SetTextColor(ROOT.kBlue)
		
		htot.SetMinimum(10e-6)
		htot.SetMaximum(10e5)
		
		ps.c.Update()   
# 		ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
# 		
# 		htot.GetXaxis().SetRangeUser(600,1200)
# 		
# 		ps.save('mass%i_%i_ZOOM' % (600, 1200), pdf=True, pdf_log=True)
# 		ps.c.Update()
		
		canvas.cd()
		pad2 = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.25)
		pad2.SetTopMargin(0)
		pad2.SetBottomMargin(0.2)
		pad2.SetGridx()
		pad2.Draw()
		pad2.cd()
		
		xax = htot.GetXaxis()
		if 'LeptonPt' != variable:
			hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';p_{T} (GeV);(fit-MC)/MC', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
		else:
			hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(NoScale)/ NoScale', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
		
		for h in [hres]:
			h.SetMarkerStyle(2)
			h.SetStats(0)
			
		for i in xrange(1, hres.GetNbinsX()+1):
			xlo = xax.GetBinLowEdge(i)
			xhi = xax.GetBinLowEdge(i+1)
			if xlo >= lo and xhi <= hi:#and xhi < 500:
				if 'LeptonPt' == variable:
					res = htot_NoScale.GetBinContent(i) / htot.GetBinContent(i)
				else:
					res = fcn.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
												
				if htot.GetBinContent(i) > 0:
					if 'LeptonPt' == variable:
						noscale = htot_NoScale.GetBinContent(i)
						noscale_err = htot_NoScale.GetBinError(i)
						scale = htot.GetBinContent(i)	
						scale_err = htot.GetBinError(i)
						err = sqrt((noscale_err/scale)*(noscale_err/scale) + (scale_err*res/scale)*(scale_err*res/scale))
						hres.SetBinContent(i, res)
						hres.SetBinError(i, err)
					else:
						hres.SetBinContent(i, res/htot.GetBinContent(i))
						hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))
		if 'LeptonPt' == variable:
			hres.SetMinimum(0)
			hres.SetMaximum(2)
		else:
			hres.SetMinimum(-1)
			hres.SetMaximum(1)
		hres.GetXaxis().SetRangeUser(lo, hi)
		hres.SetLineColor(1)
		hres.SetLineWidth(1)
		hres.SetTitle(" ")
		hres.GetYaxis().SetTitle("NoScale/Scale")
		hres.GetYaxis().SetTitleSize(20)
		hres.GetYaxis().SetTitleFont(43)
		hres.GetYaxis().SetLabelSize(0.07)
		hres.GetXaxis().SetLabelSize(0.1)
		hres.SetStats(0)
		hres.GetXaxis().SetTitleSize(20)
		if variable == 'LeptonPt':
			hres.GetXaxis().SetTitle("p_{T} [GeV]")
			l1 = ROOT.TLine(lo, 1., hi,  1.)
		else:
			hres.GetXaxis().SetTitle("M_{#mu^{+}#mu^{-}} [GeV]")
			l1 = ROOT.TLine(lo, 0., hi,  0.)
		hres.Draw('e')
		l1.Draw()
		l1.SetLineColor(ROOT.kGreen)
		
		canvas.Update()
		canvas.Print('./plots/SCALE_' + bia + '/' + variable + '/Mass_%i_%i.png' % (lo, hi))
		canvas.Print('./plots/SCALE_' + bia + '/' + variable + '/Mass_%i_%i.pdf' % (lo, hi))

		if 'LeptonPt' == variable:
			canvas.Update()			
			htot.GetXaxis().SetRangeUser(600,1200)
			htot_NoScale.GetXaxis().SetRangeUser(600,1200)
			hres.GetXaxis().SetRangeUser(600,1200)	
			canvas.Update()	
			canvas.Print('./plots/SCALE_' + bia + '/' + variable + '/Mass_%i_%i_ZOOM.png' % (lo, hi))
			canvas.Print('./plots/SCALE_' + bia + '/' + variable + '/Mass_%i_%i_ZOOM.pdf' % (lo, hi))
		
# 		ps.save('res_%i_%i' % (lo, hi), pdf=True, log=False)
# 		ps.c.Update()
# 		ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
# 		
# 		htot.GetXaxis().SetRangeUser(600,1200)
# 		
# 		ps.save('mass%i_%i_ZOOM' % (600, 1200), pdf=True, pdf_log=True)
# 		ps.c.Update()
		





		
		canvas_dati = ROOT.TCanvas(variable + "_dati", variable + "_dati", 300, 60, 1000, 650)
		pad11 = ROOT.TPad("pad1", "pad1", 0, 0.23, 1, 1.0)
		pad11.SetBottomMargin(0)
		pad11.SetGridx()
		pad11.SetLogy()
		pad11.Draw()
		pad11.cd() 
		
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
# 		if 'LeptonPt' != variable:
# 			htot.Fit(fcn, 'SREMV')
		data.Draw("same")
		legend = ROOT.TLegend(0.3, 0.75, 0.6, 0.85)
		legend.AddEntry(data, "Data", "l")
		legend.AddEntry(htot, "MC", "l")
		if 'LeptonPt' != variable:
			legend.AddEntry(fcn, "Fit MC", "l")
		bias_label.Draw()
		legend.Draw()
		htot.GetXaxis().SetRangeUser(lo, high_for_data)
		data.GetXaxis().SetRangeUser(lo, high_for_data)
		data.SetMaximum(10e5)
		data.SetMinimum(10e-6)
		ps.c.Update() 
		
# 		ps.save('data_%i_%i' % (lo, hi), pdf_log=True)
# 		
# 		data.SetMinimum(10)
# 		data.GetXaxis().SetRangeUser(600,1200)
# 		ps.save('data%i_%i_ZOOM' % (600, 1200), pdf_log=True)

		canvas_dati.cd()
		pad22 = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.23)
		pad22.SetTopMargin(0)
		pad22.SetBottomMargin(0.2)
		pad22.SetGridx()
		pad22.Draw()
		pad22.cd()
		
		data.GetXaxis().SetRangeUser(lo, high_for_data)
		if 'LeptonPt' == variable:
			hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';p_{T} (GeV);(MC-data)/data', data.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(data.GetNbinsX()+1))
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
					res = data.GetBinContent(i) / htot.GetBinContent(i)#fcn.Integral(xlo, xhi)/(xhi-xlo) - data.GetBinContent(i)
				else:
					res = data.GetBinContent(i) / htot.GetBinContent(i)
# 				print " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %f" %res
				if data.GetBinContent(i) > 0:
					if 'LeptonPt' != variable:
						data_con = data.GetBinContent(i)
						data_con_err = data.GetBinError(i)
						scale = htot.GetBinContent(i)	
						scale_err = htot.GetBinError(i)
						err = sqrt((data_con_err/scale)*(data_con_err/scale) + (scale_err*res/scale)*(scale_err*res/scale))
						hres.SetBinContent(i, res)
						hres.SetBinError(i, err)
					else:
						hres.SetBinContent(i, res/data.GetBinContent(i))
						hres.SetBinError(i, data.GetBinError(i)/data.GetBinContent(i))
				if data.GetBinContent(i) == 0:
					hres.SetBinContent(i, 1)
					hres.SetBinError(i, 0)
		if 'LeptonPt' == variable:
			hres.SetMinimum(0)
			hres.SetMaximum(2)
		else:
			hres.SetMinimum(0)
			hres.SetMaximum(2)

		hres.GetXaxis().SetRangeUser(lo, hi)
		hres.SetLineColor(1)
		hres.SetLineWidth(1)
		hres.SetTitle(" ")
		hres.GetYaxis().SetTitle("Data/MC")
		hres.GetYaxis().SetTitleSize(20)
		hres.GetYaxis().SetTitleFont(43)
		hres.GetYaxis().SetLabelSize(0.07)
		hres.GetXaxis().SetLabelSize(0.1)
		hres.SetStats(0)
		hres.GetXaxis().SetTitleSize(20)
		hres.GetXaxis().SetTitle("M_{#mu^{+}#mu^{-}} [GeV]")
		if variable == 'LeptonPt':
			hres.GetXaxis().SetTitle("p_{T} [GeV]")
			l2 = ROOT.TLine(lo, 1., hi,  1.)
		else:
			hres.GetXaxis().SetTitle("M_{#mu^{+}#mu^{-}} [GeV]")
			l2= ROOT.TLine(lo, 0., hi,  0.)
		hres.GetXaxis().SetRangeUser(lo, high_for_data)
# 		hres = hres.Rebin(12, "ciola", bin);
		ps.c.Update() 
		hres.Draw('e')
		l2.SetLineColor(ROOT.kGreen)
		l2.Draw()

		canvas_dati.Update()
		canvas_dati.Print('./plots/SCALE_' + bia + '/' + variable + '/Dati_%i_%i.png' % (lo, hi))
		canvas_dati.Print('./plots/SCALE_' + bia + '/' + variable + '/Dati_%i_%i.pdf' % (lo, hi))

		if 'LeptonPt' == variable:
			canvas.Update()			
			htot.GetXaxis().SetRangeUser(600,1200)
			data.SetMaximum(1000)
			data.SetMinimum(0.1)
			data.GetXaxis().SetRangeUser(600,1200)
			hres.GetXaxis().SetRangeUser(600,1200)	
			canvas_dati.Update()	
			canvas_dati.Print('./plots/SCALE_' + bia + '/' + variable + '/Dati_%i_%i_ZOOM.png' % (lo, hi))
			canvas_dati.Print('./plots/SCALE_' + bia + '/' + variable + '/Dati_%i_%i_ZOOM.pdf' % (lo, hi))

		
# 		ps.save('res_%i_%i_data' % (lo, hi), pdf=True, log=False)
# 		ps.c.Update() 
		
    
#    print fcn.GetProb(), fcn.GetParameter(2), fcn.GetParError(2)
#     c = r.GetCovarianceMatrix()
#     r.Print("V")
#     print " ----------------------------------------------------------------------------------====================== FINE "
#     print c

	l = [150] 
	for lo in l:
		fit_it(lo, high)
		print " ------------------------------------------------------------------------- %s" %bia

print " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
print variable, high, Zpeak
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
# draw_overlay(fsets)
