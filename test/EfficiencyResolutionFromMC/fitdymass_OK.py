#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)


ROOT.gStyle.SetOptStat(10)

# ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptFit(1111)
# ROOT.gStyle.SetOptFit(1011)

variable = 'DimuonMassVertexConstrained'
# variable = 'DimuonMassVertexConstrained_bb'
# variable = 'DimuonMassVertexConstrained_be'
# variable = 'DileptonMass_bb'
# variable = 'DileptonMass_be'
# variable = 'DileptonMass'

ps = plot_saver('plots/Fit/' + variable)

# variable = 'DileptonPt'
# variable = 'LeptonPt'
# variable = 'DileptonMass_bb'



low = fitlow =150
high = fithigh = 5000
high_for_data = 5000


int_lumi_OLD = 36238.734
int_lumi = 39484.0

# if '_bb' in variable:
# 	rescale_factor = 0.9880
# elif '_be' in variable:
# 	rescale_factor = 0.9625
# else:
# 	rescale_factor = 0.9714

rebin = 40
use_non_dy = True

masses_OLD = ["dyInclusive50", 
	 		"qcd80to120", "qcd120to170", "qcd170to300", "qcd300to470", "qcd470to600", "qcd600to800", "qcd800to1000", "qcd1000to1400", "qcd1400to1800", "qcd1800to2400", "qcd2400to3200", "qcd3200", 
	 		"Wantitop", "tW", 
	 		"ttbar_lep50to500", "ttbar_lep_500to800", "ttbar_lep_800to1200", "ttbar_lep_1200to1800", "ttbar_lep1800toInf", 
	 		"Wjets", 
	 		"WWinclusive", "WW200to600", "WW600to1200", "WW1200to2500", "WW2500", 
	 		"ZZ", "WZ", "ZZ_ext", "WZ_ext", 
	 		"dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000"]

nevents_OLD = [19385554,  6986740, 6708572, 6958708, 4150588, 3959986, 3896412, 3992112, 2999069, 396409, 397660, 399226, 391735, 
				6933094, 6952830,
				79092400, 200000, 199800, 200000, 40829, 
				29705748,
				1999000, 200000, 200000, 200000, 38969, 
				990064, 1000000, 998034, 2995828,
				100000, 100000, 100000, 98000, 96613, 100000, 100000, 100000, 100000
				]

sigmas_OLD  = [6025.2, 2762530, 471100, 117276, 7823, 648.2, 186.9, 32.3992112, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445,
			35.6, 35.6,
			87.31, 0.32611, 0.03265, 0.00305, 0.00017, 
			61526.7,
			12.178, 1.385, 0.0566, 0.0035, 0.00005,
			16.523, 47.13, 16.523, 47.13,
			1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7
			]


masses = [	"Wantitop", "tW", 
	 		"ttbar", 
	 		"WW",
	 		"ZZ", "WZ",
	 		"dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000"
	 		]

nevents = [6620324,  6723341, 
				33844772,
				1000000,
				992884, 994554,
				100000, 100000, 100000, 98000, 96613, 100000, 100000, 100000, 100000
				]

sigmas  = [35.6, 35.6,
			831.76,
			118.7,
			16.523, 47,13,
			1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7
			]


weights_OLD = [int_lumi_OLD / nev * sig for nev,sig in zip(nevents_OLD, sigmas_OLD)]
weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]

hists_OLD = []
hists_dir_OLD = '/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/'
print " ---- 2016 ----- "
for m,w in zip(masses_OLD, weights_OLD):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_OLD + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_OLD.append(h)
    
print " ---- Passo al 2017 ----- "

hists = []
hists_dir = '../DataMCSpectraComparison/mc/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists.append(h)

  
  
htot_OLD = hists_OLD[0].Clone()
for j in xrange(1, len(hists_OLD)):
    htot_OLD.Add(hists_OLD[j])
    
htot = hists[0].Clone()
for j in xrange(1, len(hists)):    
    htot.Add(hists[j])


htot_OLD.SetTitle('')
htot_OLD.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
# htot_OLD.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.

htot.SetTitle('')
htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
# htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.


htot_OLD.SetLineColor(ROOT.kBlue)
htot.SetLineColor(ROOT.kRed)


htot.Draw()
htot_OLD.Draw()

ratio = htot_OLD.Clone()
htot_for_ratio = htot.Clone()

def fit_it(lo, hi):

    
    print " ----------------------------------------------------------------------------------====================== INIZIO "
    Ymin = 10e-5
    Ymax = 10e6
	   
    fcn = ROOT.TF1("fcn", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn.SetLineColor(ROOT.kBlue)
    

    fcn_systematic = ROOT.TF1("fcn_systematic", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn_systematic.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn_systematic.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn_systematic.SetLineColor(ROOT.kRed)
    
    
    htot_OLD.Fit(fcn, 'SREMV same')
    htot_OLD.SetMinimum(Ymin)
    htot_OLD.SetMaximum(Ymax)
    ps.c.Update()

    htot.Fit(fcn_systematic, 'SREMV same')
    htot.SetMinimum(Ymin)
    htot.SetMaximum(Ymax)
    ps.c.Update()


    ss = htot.GetListOfFunctions().FindObject("stats")
    ss.SetName("Const")
    ss.SetX1NDC(0.7)
    ss.SetY1NDC(0.32)
    ss.SetY2NDC(0.62)
#     ss.SetOptStat(10)
    ss.SetTextColor(ROOT.kRed)
#     ss.SetOptFit(1111)
    ps.c.Update()
    
    s = htot_OLD.GetListOfFunctions().FindObject("stats")
    s.SetName("Const")
    s.SetX1NDC(0.7)
    s.SetY1NDC(0.65)
    s.SetY2NDC(0.95)
#     s.SetOptStat(10)
    s.SetTextColor(ROOT.kBlue)
#     s.SetOptFit(1111)
    ps.c.Update()

    s.Draw()
    ss.Draw()
    
    legend = ROOT.TLegend(0.3, 0.7, 0.6, 0.85)
    legend.AddEntry(htot_OLD, "2016", "l")
    legend.AddEntry(htot, "2017", "l")
    legend.Draw()

    ps.c.Update()
    
    ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
        
    
    ratio.Divide(htot_for_ratio)
#     ratio = ratio - 1
    ratio.SetLineColor(1)
    ratio.SetLineWidth(1)
#     for i in xrange(1, ratio.GetNbinsX()+1):
#     	bin_old = ratio.GetBinContent(i)
#     	ratio.SetBinContent(i, bin_old -1)

#     ratio.GetYaxis().SetTitle("Old / New - 1")
    ratio.GetYaxis().SetTitle("Old / New")
    ratio.SetStats(1)
#     ratio.SetOptStats(111)
    line = ROOT.TLine(150, 0, 5000, 0)
    line.SetLineColor(2)
    line.SetLineWidth(1)
    ratio.Draw()
#     ratio.GetXaxis().SetRangeUser(0,240)
    ps.c.Update()
    print "========================================================================================"    
    fcn_ratio = ROOT.TF1("fcn_ratio", "[0] + [1]*x")
#     fcn_ratio = ROOT.TF1("fcn_ratio", "[0] + [1] * x + [2] * x*x")
    print "ok"
    ratio.Fit(fcn_ratio)#, 'SREMV same')
    sss = ratio.GetListOfFunctions().FindObject("stats")
    sss.SetX1NDC(0.73)
    sss.SetY1NDC(0.69)
    sss.SetY2NDC(0.94)
    sss.SetOptStat(10)
    sss.Draw("same")
    ps.c.Update()
#     ratio.GetYaxis().SetRangeUser(-1, 1)
    ratio.GetYaxis().SetRangeUser(0, 2)
    line.Draw()

      
    ps.save('res_%i_%i' % (lo, hi), pdf=True, log=False)
#     
    ps.c.Update()
#        
#        

	

#l = range(200, 2200, 200)
# l = [60, 120, 140, 160, 200]
l = [70]
for lo in l:
    print " ----------------------------------------------------------------------------------->>>>>>>>>>>>>>>>>>>>> %d " %lo
    fit_it(lo, high)
    
print " ------------------------------------------------------------------------- Passo al Draw "
print variable, high,int_lumi, use_non_dy, rebin#, bias
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
# draw_overlay(fsets)
