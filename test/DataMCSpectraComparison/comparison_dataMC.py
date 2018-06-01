import sys
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import overall_prescale


ROOT.TH1.AddDirectory(False)



histos = {}
histos_7 = {}
mumu_MC = '/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/ana_datamc_%s.root'
mumu_Dati = '/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data_OK/NoScale_YesEtaCut_Run2016MuonsOnly/ana_datamc_data.root'
# mumu_Dati = "/afs/cern.ch/work/f/ferrico/private/ZPrime_2017/CMSSW_9_2_6/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data_OK/Run2017MuonsOnly/ana_datamc_data.root"

mumu_scale = 1
# lumi = 14340.494
lumi = 36235.493 #36295.39 #37197 36295 #
# histogram = ['DileptonMass', 'DimuonMassVertexConstrained']#'DileptonMass_bb', 'DileptonPt']
# histogram = ['LeptonPt']
histogram = ['DimuonMassVertexConstrained_bb', 'DimuonMassVertexConstrained_be']#, 'DimuonMassVertexConstrained_be']
#,  'DileptonMass_be',  'DimuonMassVertexConstrained_be', 'LeptonPt']
# 
Xmax = 3000
Xmin = 60

tt_binned = True;

if tt_binned:
	tt_sample = 'Binned'
else:
	tt_sample = 'NoBinned'

rebin = 20
for mumu_histogram in histogram:
	if 'Pt' in mumu_histogram:
		rebin = 20

	if '_bb' in mumu_histogram:
		Zpeak = 0.9638
	elif '_be' in mumu_histogram:
		Zpeak = 0.9688
	else:
		Zpeak = 0.9610
		
	MC_Stack = ROOT.THStack('MC', '')
	MC_Stack_cum = ROOT.THStack('MC_cum', '')
	dati = ROOT.TH1F()
	h = ROOT.TH1F()
	h_cum = ROOT.TH1F()

	for sample in samples:
	

 		file_dati = ROOT.TFile(mumu_Dati)
		dati = file_dati.Our2016MuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
		dati.Rebin(rebin)
		dati.SetMarkerColor(1)
		dati.SetLineColor(1)
		dati.SetLineWidth(2)
		
		f = ROOT.TFile(mumu_MC % sample.name)
		
		h = f.Our2016MuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
		h.Scale(sample.partial_weight * lumi * Zpeak )
		print sample.name, sample.partial_weight, sample.partial_weight * lumi * Zpeak
		h.Rebin(rebin)
		h_cum = f.Our2016MuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
		h_cum.Scale(sample.partial_weight * lumi * Zpeak )
		h_cum.Rebin(rebin)		

		nbins = h_cum.GetNbinsX()
		contenuto_mc = 0
		for i in xrange(0, nbins):
# 			contenuto_mc += h.GetBinContent(nbins - i)
			contenuto_mc = h.Integral(i, nbins)
			h_cum.SetBinContent(i, contenuto_mc)
# 			if 'dy120to200' in sample.name or 'dy200to400' in sample.name:
# 				if i*h.GetBinWidth(i) < 300:
# 					print sample.name, i*h.GetBinWidth(i), contenuto_mc, h.Integral(i, nbins) - h.Integral(i-1, nbins), h_cum.GetBinContent(i)
				

# 		h.SetLineWidth(1)

		if 'dy' in sample.name:
			if 'Inclusive' in sample.name:
				h.SetLineColor(ROOT.kYellow)
				h.SetFillColor(ROOT.kYellow)
				h.SetMarkerColor(ROOT.kYellow)
				h_cum.SetLineColor(ROOT.kYellow)
				h_cum.SetFillColor(ROOT.kYellow)
				h_cum.SetMarkerColor(ROOT.kYellow)
			else:
				h.SetLineColor(3)
				h.SetFillColor(3)
				h.SetMarkerColor(3)
				h_cum.SetLineColor(3)
				h_cum.SetFillColor(3)
				h_cum.SetMarkerColor(3)
				for i in range(0, h.GetNbinsX()):
					width = h.GetBinWidth(i)
					if(i*width > 60 and i*width < 120):
						old = h.GetBinContent(i)
						new = old * 0.98
						h.SetBinContent(i, new)

						old = h_cum.GetBinContent(i)
						new = old * 0.98
						h_cum.SetBinContent(i, new)
					if(i*width > 120 and i*width < 150):
						old = h.GetBinContent(i)
						new = old * 1.24
						h.SetBinContent(i, new)

						old = h_cum.GetBinContent(i)
						new = old * 1.24
						h_cum.SetBinContent(i, new)
					
		elif 'ZZ' in sample.name:
			h.SetFillColor(2)
			h.SetLineColor(2)
			h.SetMarkerSize(2)
			h.SetMarkerColor(2)
			h_cum.SetFillColor(2)
			h_cum.SetLineColor(2)
			h_cum.SetMarkerSize(2)
			h_cum.SetMarkerColor(2)

		elif 'WZ' in sample.name:
			h.SetFillColor(ROOT.kRed+3)
			h.SetLineColor(ROOT.kRed+3)
			h.SetMarkerSize(ROOT.kRed+3)
			h.SetMarkerColor(ROOT.kRed+3)
			h_cum.SetFillColor(ROOT.kRed+3)
			h_cum.SetLineColor(ROOT.kRed+3)
			h_cum.SetMarkerSize(ROOT.kRed+3)
			h_cum.SetMarkerColor(ROOT.kRed+3)
			
		elif 'WW' in sample.name:
			h.SetFillColor(ROOT.kRed-10)
			h.SetLineColor(ROOT.kRed-10)
			h.SetMarkerSize(ROOT.kRed-10)
			h.SetMarkerColor(ROOT.kRed-10)
			h_cum.SetFillColor(ROOT.kRed-10)
			h_cum.SetLineColor(ROOT.kRed-10)
			h_cum.SetMarkerSize(ROOT.kRed-10)
			h_cum.SetMarkerColor(ROOT.kRed-10)
			
# 			 or 'WW' in sample.name or sample.name == 'WZ':
		elif sample.name == 'Wantitop' or sample.name == 'tW':
			h.SetFillColor(ROOT.kBlue+2)
			h.SetLineColor(ROOT.kBlue+2)
			h.SetMarkerSize(ROOT.kBlue+2)
			h.SetMarkerColor(ROOT.kBlue+2)
			h_cum.SetFillColor(ROOT.kBlue+2)
			h_cum.SetLineColor(ROOT.kBlue+2)
			h_cum.SetMarkerSize(ROOT.kBlue+2)
			h_cum.SetMarkerColor(ROOT.kBlue+2)
			
		elif 'qcd' in sample.name or 'jet' in sample.name:
			h.SetLineColor(7)
			h.SetMarkerSize(7)
			h.SetMarkerColor(7)
			h.SetFillColor(7)
			h_cum.SetLineColor(7)
			h_cum.SetMarkerSize(7)
			h_cum.SetMarkerColor(7)
			h_cum.SetFillColor(7)
			
		elif 'tt' in sample.name:
			h.SetFillColor(4)
			h.SetLineColor(4)
			h.SetMarkerSize(4)
			h.SetMarkerColor(4)
			h_cum.SetFillColor(4)
			h_cum.SetLineColor(4)
			h_cum.SetMarkerSize(4)
			h_cum.SetMarkerColor(4)


# 		h.GetYaxis().SetRangeUser(10e-4, 2*10e4)
# 		h_cum.GetYaxis().SetRangeUser(10e-4, 2*10e4)
		MC_Stack.Add(h, "HIST")
		MC_Stack_cum.Add(h_cum, "HIST")

		print sample.name, h.Integral(120,400)
    		
 
 	totale_80X = MC_Stack.GetStack().Last()
 	totale_80X_cum = MC_Stack_cum.GetStack().Last()
# 	totale_80X.SetLineColor(2)
	
	canvas = ROOT.TCanvas(mumu_histogram, mumu_histogram, 210,45,1050,750)
	canvas.SetGrid()
	pad1 = ROOT.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
	pad1.SetBottomMargin(0)
	pad1.SetGridx()
	pad1.Draw()
	pad1.cd() 

	pad1.SetLogy()
	
	label_histogram = ROOT.TPaveLabel(0.30, 0.775, 0.60, 0.875, mumu_histogram, 'brNDC')
	label_histogram.SetTextFont(42)
	label_histogram.SetTextSize(0.5)
	label_histogram.SetBorderSize(0)
	label_histogram.SetFillColor(0)
	label_histogram.SetFillStyle(0)
	
	label_luminosity = ROOT.TPaveLabel(0.4, 0.675, 0.50, 0.775, "\int L\, dt = 36.238 fb^{-1}",'brNDC')
	label_luminosity.SetTextFont(42)
	label_luminosity.SetTextSize(0.5)
	label_luminosity.SetBorderSize(0)
	label_luminosity.SetFillColor(0)
	label_luminosity.SetFillStyle(0)

	
	#totale_80X.Scale(1/totale_80X.Integral())
	#totale_80X.Draw()
	MC_Stack.Draw()
	dati.Draw("same")
# 	totale_80X.Draw("same")
	
	#totale_80X.Rebin(rebin)
	#h.Rebin(rebin)
	#totale_80X.GetXaxis().SetTitle("Mass [GeV/c^{2}]")
	titolo_yAxis = "event / %s GeV" % rebin
	totale_80X.GetYaxis().SetTitle(titolo_yAxis)
	totale_80X.SetStats(0)
# 	h.SetStats(0)
# 	h_cum.SetStats(0)
	dati.SetStats(0)
	

# 	Xmax = totale_80X.GetXaxis().GetXmax()
# 	Xmin = totale_80X.GetXaxis().GetXmin()

 	if 'LeptonPt' == mumu_histogram:
 		Xmin = 0
 		
 	if 'LeptonPt' == mumu_histogram:
 		Xmax = 1900
	
	#MC_Stack.GetXaxis().SetRangeUser(Xmin, Xmax)
# 	h.GetYaxis().SetRangeUser(10e-4, 2*10e4)
 	MC_Stack.GetXaxis().SetRangeUser(Xmin, Xmax)

 	if 'ass' in mumu_histogram:
	 	MC_Stack.GetXaxis().SetTitle("m_#mu#mu [GeV]")

 		
#   	MC_Stack.GetYaxis().SetRangeUser(10e-4, 2*10e4)
#  	dati.GetXaxis().SetRangeUser(Xmin, Xmax)
#   	dati.GetYaxis().SetRangeUser(10e-4, 2*10e4)
  	MC_Stack.SetMinimum(0.000001)
  	MC_Stack.SetMaximum(500000)

#   	dati.GetYaxis().SetMinimum(10e-4)
#   	dati.GetYaxis().SetMaximum(2*10e4)

  	
# 	legend = ROOT.TLegend(0.7,0.8,0.9,0.9)
# 	legend.AddEntry(MC_Stack, "MC", "f")
# 	legend.AddEntry(totale_80X,"MC", "lep")	
# 	legend.AddEntry(dati,"data: %.3f fb-1" % lumi/1000, "lep")	
# 	legend.Draw()

# 	label_histogram.Draw()
# 	label_luminosity.Draw()
	
	canvas.cd()
	pad2 = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.25)
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.2)
	pad2.SetGridx()
	pad2.Draw()
	pad2.cd()  
	nome_ratio = '%s_ratio' % mumu_histogram
	#nBin = totale_80X.GetNbinsX()#/rebin
	ratio = dati.Clone()
	ratio.Divide(totale_80X)
	ratio.SetLineColor(1)
	ratio.SetLineWidth(1)
	ratio.SetTitle(" ")
	ratio.GetYaxis().SetTitle("Data/MC")
	ratio.GetYaxis().SetTitleSize(20)
	ratio.GetYaxis().SetTitleFont(43)
	ratio.GetYaxis().SetLabelSize(0.07)
	ratio.GetXaxis().SetLabelSize(0.1)
	ratio.SetStats(0)
	ratio.GetXaxis().SetRangeUser(Xmin, Xmax)
	ratio.GetYaxis().SetRangeUser(0, 2)	
	ratio.GetXaxis().SetTitleSize(50)
	ratio.GetXaxis().SetTitle("M_{#mu^{+}#mu^{-}} [GeV]")
	line = ROOT.TLine(Xmin, 1, Xmax, 1)
	line.SetLineColor(2)
	line.SetLineWidth(1)
	ratio.Draw("ep")
	line.Draw("same")


	canvas.Update()
	canvas.Print('./Comparison_dataMC/' + tt_sample + '/%s.png' % mumu_histogram)
	canvas.Print('./Comparison_dataMC/' + tt_sample + '/%s.pdf' % mumu_histogram)
	
	pad1.SetLogx()
	pad2.SetLogx()
	canvas.Print('./Comparison_dataMC/' + tt_sample + '/%s_Log.png' % mumu_histogram)
	canvas.Print('./Comparison_dataMC/' + tt_sample + '/%s_Log.pdf' % mumu_histogram)

 	nome_file_output ='./Comparison_dataMC/' + mumu_histogram + '_.root'
	prova = ROOT.TFile(nome_file_output , 'RECREATE')
	MC_Stack.Write()
	h.Write()
	

	cumulative_res = dati.Clone()	
	canvas_cum = ROOT.TCanvas(mumu_histogram + ' cumulative', mumu_histogram + ' cumulative', 210,45,1050,750)
	canvas_cum.SetGrid()
	pad11 = ROOT.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
	pad11.SetBottomMargin(0)
	pad11.SetGridx()
	pad11.Draw()
	pad11.cd()
	pad11.SetLogy()
# 	totale_80X_cum = totale_80X.Clone()
	dati_cum = dati.Clone()
	nbins = dati_cum.GetNbinsX()
	contenuto_dati = 0
	count_mc = 0
	count_dati = 0
	for i in xrange(0, nbins):
		contenuto_dati = dati.Integral(i, nbins)
		dati_cum.SetBinContent(i, contenuto_dati)
# 		contenuto_dati += dati.GetBinContent(nbins - i)
# 		dati_cum.SetBinContent(nbins-i, contenuto_dati)
		
		
# 		for sample in samples:
# 			contenuto_mc += MC_Stack.GetStack().at(i).GetBinContent(nbins - i)
# 
# 			MC_Stack_cum.SetBinContent(nbins-i, contenuto_mc)
# 		if i * totale_80X.GetBinWidth(i) >= 1700 and i * totale_80X.GetBinWidth(i) < 2300:
# 			count_mc += totale_80X.GetBinContent(i)
# 		if i * dati.GetBinWidth(i) >= 1700 and i * dati.GetBinWidth(i) < 2300:
# 			count_dati += dati.GetBinContent(i)
# 			if 'DileptonMass_bb' == mumu_histogram or 'DimuonMassVertexConstrained_bb' == mumu_histogram:
# 				print i * dati.GetBinWidth(i), dati.GetBinContent(nbins - i)
		
				
# 		totale_80X_cum.SetBinContent(i, totale_80X.Integral(0, nbins - i))
# 		dati_cum.SetBinContent(i, dati.Integral(0, nbins - i))

	MC_Stack_cum.Draw()
	dati_cum.Draw("same")
 	MC_Stack_cum.GetXaxis().SetRangeUser(Xmin, Xmax)
  	MC_Stack_cum.SetMinimum(0.1)
  	MC_Stack_cum.SetMaximum(500000)
  	
	canvas_cum.cd()
	pad22 = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.25)
	pad22.SetTopMargin(0)
	pad22.SetBottomMargin(0.2)
	pad22.SetGridx()
	pad22.Draw()
	pad22.cd()
	cumulative_res = dati_cum.Clone()
	cumulative_res.Divide(totale_80X_cum)
	cumulative_res.SetLineColor(1)
	cumulative_res.SetLineWidth(1)
	cumulative_res.SetTitle(" ")
	cumulative_res.GetYaxis().SetTitle("Data/MC")
	cumulative_res.GetYaxis().SetTitleSize(20)
	cumulative_res.GetYaxis().SetTitleFont(43)
	cumulative_res.GetYaxis().SetLabelSize(0.07)
	cumulative_res.GetXaxis().SetLabelSize(0.1)
	cumulative_res.SetStats(0)
	cumulative_res.GetXaxis().SetRangeUser(Xmin, Xmax)
	cumulative_res.GetYaxis().SetRangeUser(0, 2)	
	cumulative_res.GetXaxis().SetTitleSize(20)
	cumulative_res.GetXaxis().SetTitle("M_{#mu^{+}#mu^{-}} [GeV]")  
	cumulative_res.Draw("ep")
	line.Draw("same")
	
	if 'DileptonMass_bb' == mumu_histogram:
		print mumu_histogram, count_mc, count_dati, 
	if 'DimuonMassVertexConstrained_bb' == mumu_histogram:
		print mumu_histogram, count_mc, count_dati, 
	canvas.Update()
	canvas_cum.Print('./Comparison_dataMC/' + tt_sample + '/%s_cumulative.png' % mumu_histogram)
	canvas_cum.Print('./Comparison_dataMC/' + tt_sample + '/%s_cumulative.pdf' % mumu_histogram)
	
	pad11.SetLogx()
	pad22.SetLogx()
	canvas_cum.Print('./Comparison_dataMC/' + tt_sample + '/%s_cumulative_Log.png' % mumu_histogram)
	canvas_cum.Print('./Comparison_dataMC/' + tt_sample + '/%s_cumulative_Log.pdf' % mumu_histogram)

	prova.Close()