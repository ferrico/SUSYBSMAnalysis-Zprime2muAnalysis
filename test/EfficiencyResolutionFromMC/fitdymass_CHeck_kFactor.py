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

ps = plot_saver('plots/CANCELLA_SUBITO/' + variable)

# variable = 'DileptonPt'
# variable = 'LeptonPt'
# variable = 'DileptonMass_bb'



low = fitlow =150
high = fithigh = 5000
high_for_data = 5000


int_lumi = 36238.734 #36295.39#

# if '_bb' in variable:
# 	rescale_factor = 0.9880
# elif '_be' in variable:
# 	rescale_factor = 0.9625
# else:
# 	rescale_factor = 0.9714

rebin = 40
use_non_dy = True

masses  = ['dy50to120', 'dy120to200', 'dy200to400', 'dy400to800', 'dy800to1400', 'dy1400to2300', 'dy2300to3500', 'dy3500to4500', 'dy4500to6000']
nevents = [2977600, 100000, 100000,  98400, 100000, 95106, 100000, 100000, 100000]
sigmas  = [  1975, 19.32, 2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7]

weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]
# weights = [rescale_factor * int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]

hists_DY = []
# hists_dir_DY = '../DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/'
hists_dir_DY = '/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_DY + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
#     print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_DY.append(h)


hists_DY_systematic = []
hists_dir_DY_systematic = '../DataMCSpectraComparison/mc/Old_KFactor/'#NEW_KFACTOR/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_DY_systematic + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
#     print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_DY_systematic.append(h)


hists_dir = '/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/' 
if use_non_dy:
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
    non_dy_samples = [          WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500,
                                                WZ,
#                                               WZ_skim,
                                                        ZZ,
#                                                       ZZ_skim,
                                                ZZ_ext,
#                                               ZZ_ext_skim,
                                                WZ_ext,
                                                Wantitop, tW,
                                                Wjets,
                                                ttbar_lep50to500,
#                                               ttbar_lep50to500_OLD, 
                                                        ttbar_lep_500to800,
                                                        ttbar_lep_800to1200,
                                                        ttbar_lep_1200to1800,
                                                        ttbar_lep1800toInf,
                                                        #                                               qcd50to80, 
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
        w = sample.partial_weight * int_lumi# * rescale_factor
        h = d.Get(variable).Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        h.Draw()
        ps.save('rawmass_%s' % sample.name)
        print sample.name, w
        hists_DY.append(h)
        hists_DY_systematic.append(h)
  
  
htot_DY = hists_DY[0].Clone()
for j in xrange(1, len(hists_DY)):
    htot_DY.Add(hists_DY[j])
    
htot_DY_systematic = hists_DY_systematic[0].Clone()
for j in xrange(1, len(hists_DY_systematic)):    
    htot_DY_systematic.Add(hists_DY_systematic[j])


htot_DY.SetTitle('')
htot_DY.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
# htot_DY.GetXaxis().SetTitle('lepton p_{T} (GeV)')
htot_DY.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.

htot_DY.SetLineColor(ROOT.kBlue)
htot_DY_systematic.SetLineColor(ROOT.kRed)



# htot_DY.SetMinimum(10e2)
# htot_DY_systematic.SetMinimum(10e2)
# htot_DY.SetMaximum(10e4)
# htot_DY.GetXaxis().SetRangeUser(150,500)


ratio = htot_DY.Clone()
htot_DY_systematic_2 = htot_DY_systematic.Clone()

def fit_it(lo, hi):

    
    print " ----------------------------------------------------------------------------------====================== INIZIO "
    Ymin = 10e-5
    Ymax = 10e5
	   
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
    
    htot_DY.Fit(fcn, 'SREMV same')
    htot_DY.SetMinimum(Ymin)
    htot_DY.SetMaximum(Ymax)
    ps.c.Update()

    htot_DY_systematic.Fit(fcn_systematic, 'SREMV same')
    htot_DY_systematic.SetMinimum(Ymin)
    htot_DY_systematic.SetMaximum(Ymax)
    ps.c.Update()


    ss = htot_DY_systematic.GetListOfFunctions().FindObject("stats")
    ss.SetName("Const")
    ss.SetX1NDC(0.7)
    ss.SetY1NDC(0.32)
    ss.SetY2NDC(0.62)
#     ss.SetOptStat(10)
    ss.SetTextColor(ROOT.kRed)
#     ss.SetOptFit(1111)
    ps.c.Update()
    
    s = htot_DY.GetListOfFunctions().FindObject("stats")
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
    legend.AddEntry(htot_DY, "MC", "l")
    legend.AddEntry(htot_DY_systematic, "MC - New kFactor", "l")
    legend.Draw()

    ps.c.Update()
    
    ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
        
    
    ratio.Divide(htot_DY_systematic_2)
#     ratio = ratio - 1
    ratio.SetLineColor(1)
    ratio.SetLineWidth(1)
    for i in xrange(1, ratio.GetNbinsX()+1):
    	bin_old = ratio.GetBinContent(i)
    	ratio.SetBinContent(i, bin_old -1)
#     ratio.SetTitle(" ")
    ratio.GetYaxis().SetTitle("Old / New - 1")
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
    ratio.GetYaxis().SetRangeUser(-1, 1)
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
