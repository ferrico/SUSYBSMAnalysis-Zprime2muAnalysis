import sys
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import overall_prescale


ROOT.TH1.AddDirectory(False)



histos = {}
histos_7 = {}
mumu_MC = '/afs/cern.ch/work/f/ferrico/private/Codice_ZPrime_8_NoTrigger/CMSSW_8_0_3_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/mc/ana_datamc_%s.root'
mumu_Dati = '/afs/cern.ch/work/f/ferrico/private/Codice_ZPrime_8/CMSSW_8_0_3_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data/Run2015MuonsOnly/ana_datamc_data.root'
mumu_scale = 1
rebin = 10
lumi = 6283.199#4079.344
histogram = ['DimuonMassVertexConstrained', 'DimuonMassVertexConstrained_bb', 'DimuonMassVertexConstrained_ne', 'DimuonMassVertexConstrained_pe']
for mumu_histogram in histogram:
	
	MC_Stack = ROOT.THStack('MC', '')
	dati = ROOT.TH1F()
	h = ROOT.TH1F()

	for sample in samples:
	

 		file_dati = ROOT.TFile(mumu_Dati)
		dati = file_dati.Our2012MuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
		dati.Rebin(rebin)
		dati.SetMarkerColor(1)
		dati.SetLineColor(1)
		dati.SetLineWidth(2)


		
		f = ROOT.TFile(mumu_MC % sample.name)
		h = f.Our2012MuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
		h.Scale(sample.partial_weight * lumi)
		h.Rebin(rebin)
		h.SetLineWidth(1)

		if 'dy' in sample.name:
			#h.SetFillColor(7)
			h.SetLineColor(7)
			h.SetMarkerColor(7)
			
		elif sample.name == 'ZZ' or 'WW' in sample.name or sample.name == 'WZ':
			#h.SetFillColor(ROOT.kBlue+2)
			h.SetLineColor(ROOT.kBlue+2)
			h.SetMarkerSize(ROOT.kBlue+2)
			h.SetMarkerColor(ROOT.kBlue+2)
		elif 'qcd' in sample.name or sample.name == 'Wantitop' or sample.name == 'tW':
			#h.SetFillColor(ROOT.kYellow+2)
			h.SetLineColor(ROOT.kYellow+2)
			h.SetMarkerSize(ROOT.kYellow+2)
			h.SetMarkerColor(ROOT.kYellow+2)
			
		MC_Stack.Add(h)
    		
 
 	totale_80X = MC_Stack.GetStack().Last()
# 	totale_80X.SetLineColor(2)
	
	canvas = ROOT.TCanvas(mumu_histogram, mumu_histogram, 210,45,1050,750)
	#canvas.SetOptStat(0)
	pad1 = ROOT.TPad("pad1", "pad1", 0, 0.23, 1, 1.0)
	pad1.SetBottomMargin(0)
	pad1.SetGridx()
	pad1.Draw()
	pad1.cd() 
		

	pad1.SetLogy()
	
	t = ROOT.TPaveLabel(0.50, 0.425, 0.90, 0.525, mumu_histogram, 'brNDC')
	t.SetTextFont(42)
	t.SetTextSize(0.5)
	t.SetBorderSize(0)
	t.SetFillColor(0)
	t.SetFillStyle(0)

	#totale_80X.Scale(1/totale_80X.Integral())
	#totale_80X.Draw()
	dati.Draw("")
	MC_Stack.Draw("same")

	#totale_80X.Rebin(rebin)
	#h.Rebin(rebin)
	#totale_80X.GetXaxis().SetTitle("Mass [GeV/c^{2}]")
	titolo_yAxis = "event / %s GeV" % rebin
	#totale_80X.GetYaxis().SetTitle(titolo_yAxis)
	#totale_80X.SetStats(0)
	h.SetStats(0)
	

# 	Xmax = totale_80X.GetXaxis().GetXmax()
# 	Xmin = totale_80X.GetXaxis().GetXmin()
	Xmin = 0
	Xmax = 3500
	
	#MC_Stack.GetXaxis().SetRangeUser(Xmin, Xmax)
	#h.GetYaxis().SetRangeUser(10e-2, 2*10e4)
 	dati.GetXaxis().SetRangeUser(Xmin, Xmax)
  	dati.GetYaxis().SetRangeUser(10e-6, 2*10e4)
	
	#totale_80X.GetXaxis().SetRangeUser(Xmin,Xmax)
	legend = ROOT.TLegend(0.7,0.8,1.05,0.9)
	legend.AddEntry(MC_Stack, "MC", "f")
	#legend.AddEntry(totale_80X,"MC", "f")	
	ciao = lumi/1000
	legend.AddEntry(dati,"data: %.3s fb-1" % ciao, "lep")	
	legend.Draw()
	t.Draw()

	canvas.cd()
	pad2 = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.23)
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
	ratio.GetYaxis().SetTitle("ratio: dati/MC")
	ratio.GetYaxis().SetTitleSize(20)
	ratio.GetYaxis().SetTitleFont(43)
	ratio.GetYaxis().SetLabelSize(0.07)
	ratio.GetXaxis().SetLabelSize(0.1)
	ratio.SetStats(0)
	ratio.GetXaxis().SetRangeUser(Xmin, Xmax)
	ratio.GetYaxis().SetRangeUser(0, 2)	
	line = ROOT.TLine(Xmin, 1, Xmax, 1)
	line.SetLineColor(2)
	line.SetLineWidth(1)
	ratio.Draw("ep")
	line.Draw("same")


	canvas.Update()
	canvas.Print('./Comparison_dataMC/%s.png' % mumu_histogram)
 	nome_file_output ='./Comparison_dataMC/' + mumu_histogram + '_.root'
	prova = ROOT.TFile(nome_file_output , 'RECREATE')
	MC_Stack.Write()
	h.Write()
	prova.Close()
