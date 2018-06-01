#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
# ROOT.gStyle.SetOptFit(111)

# variable = 'DimuonMassVertexConstrained'
# variable = 'DimuonMassVertexConstrained_bb'
category = 'Barrel - Barrel'

# variable = 'DileptonPt'
# variable = 'LeptonPt'
variable = 'DileptonMass_bb'


# category = 'BE + EE'
# variable = 'DimuonMassVertexConstrained_be'

low = fitlow = 150
high = fithigh = 2500
high_for_data = 2500

# ps = plot_saver('plots/fitdymass/split'+ variable)
# bias = '0p15'
# bias = '0p1'
# bias = '0p05'
# bias = ''


# ps = plot_saver('plots/SCALE_'+ bias + '/' + variable)
ps = plot_saver('plots/SCALE_sovrapposizione' + variable)

int_lumi = 36238.734

if '_bb' in variable:
	rescale_factor = 0.9880
elif '_be' in variable:
	rescale_factor = 0.9625
else:
	rescale_factor = 0.9714
rebin = 40
use_non_dy = True

# Masses = [50, 120, 200, 400, 800, 1400, 2300, 3500, 4500]


masses  = ['dy50to120', 'dy120to200', 'dy200to400', 'dy400to800', 'dy800to1400', 'dy1400to2300', 'dy2300to3500', 'dy3500to4500', 'dy4500to6000']

nevents = [2977600, 100000, 100000,  98400, 100000, 95106, 100000, 100000, 100000]

sigmas  = [  1975, 19.32, 2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7]

weights = [rescale_factor * int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]
#weights = [x/weights[-1] for x in weights]

# hists_dir_DY = '../DataMCSpectraComparison/mc/DY_scale_' + bias + '/'

hists_DY_0p01 = []
hists_dir_DY_0p01 = '../DataMCSpectraComparison/mc/DY_kFactor_Scale_0p01/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_DY_0p01 + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_DY_0p01.append(h)

hists_DY_0p05 = []
hists_dir_DY_0p05 = '../DataMCSpectraComparison/mc/DY_kFactor_Scale_0p05/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_DY_0p05 + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_DY_0p05.append(h)

hists_DY_0p1 = []
hists_dir_DY_0p1 = '../DataMCSpectraComparison/mc/DY_kFactor_Scale_0p1/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_DY_0p1 + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_DY_0p1.append(h)

hists_DY_0p15 = []
hists_dir_DY_0p15 = '../DataMCSpectraComparison/mc/DY_kFactor_Scale_0p15/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_DY_0p15 + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_DY_0p15.append(h)

hists_DY_0p2 = []
hists_dir_DY_0p2 = '../DataMCSpectraComparison/mc/DY_kFactor_Scale_0p2/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_DY_0p2 + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists_DY_0p2.append(h)

hists_dir = '../DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/' 
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
        ps.save('rawmass_%s' % sample.name)
        print sample.name, w
        hists_DY_0p01.append(h)
        hists_DY_0p05.append(h)
        hists_DY_0p1.append(h)
        hists_DY_0p15.append(h)
        hists_DY_0p2.append(h)
  
  
htot_DY_0p01 = hists_DY_0p01[0].Clone()
htot_DY_0p05 = hists_DY_0p05[0].Clone()
htot_DY_0p1 = hists_DY_0p1[0].Clone()
htot_DY_0p15 = hists_DY_0p15[0].Clone()
htot_DY_0p2 = hists_DY_0p2[0].Clone()

for j in xrange(1, len(hists_DY_0p01)):
    htot_DY_0p01.Add(hists_DY_0p01[j])
    htot_DY_0p05.Add(hists_DY_0p05[j])
    htot_DY_0p1.Add(hists_DY_0p1[j])
    htot_DY_0p15.Add(hists_DY_0p15[j])
    htot_DY_0p2.Add(hists_DY_0p2[j])

htot_DY_0p01.SetTitle('')
# htot_DY_0p01.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
htot_DY_0p01.GetXaxis().SetTitle('lepton p_{T} (GeV)')
htot_DY_0p01.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.

htot_DY_0p01.SetLineColor(ROOT.kBlue)
htot_DY_0p05.SetLineColor(ROOT.kRed)
htot_DY_0p1.SetLineColor(ROOT.kGreen)
htot_DY_0p15.SetLineColor(ROOT.kBlack)
htot_DY_0p2.SetLineColor(ROOT.kMagenta)
htot_DY_0p01.Draw()
htot_DY_0p05.Draw("same")
htot_DY_0p1.Draw("same")
htot_DY_0p15.Draw("same")
htot_DY_0p2.Draw("same")
legend = ROOT.TLegend(0.6, 0.75, 0.9, 0.95)
legend.AddEntry(htot_DY_0p01, "MC: DY with 0.01 bias", "l")
legend.AddEntry(htot_DY_0p05, "MC: DY with 0.05 bias", "l")
legend.AddEntry(htot_DY_0p1, "MC: DY with 0.1 bias", "l")
legend.AddEntry(htot_DY_0p15, "MC: DY with 0.15 bias", "l")
legend.AddEntry(htot_DY_0p2, "MC: DY with 0.2 bias", "l")
legend.Draw()

def fit_it(lo, hi):
    
    print " ----------------------------------------------------------------------------------====================== INIZIO "
    
#     htot_DY_0p01.Draw()
#     ps.c.Update()
#     htot_DY_0p05.Draw("same")
#     ps.c.Update()
#     htot_DY_0p1.Draw("same")
#     ps.c.Update()
#     htot_DY_0p15.Draw("same")
#     ps.c.Update()
#     htot_DY_0p2.Draw("same")
#     ps.c.Update()

    fcn_0p01 = ROOT.TF1("fcn_0p01", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn_0p01.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn_0p01.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn_0p01.SetLineColor(ROOT.kBlue)
    htot_DY_0p01.Fit(fcn_0p01, 'SREMV same')
#     s = htot_DY_0p01.GetListOfFunctions().FindObject("stats")
#     s.SetName("Const")
#     s.SetX1NDC(0.73)
#     s.SetY1NDC(0.69)
#     s.SetY2NDC(0.94)
#     s.SetOptStat(10)
#     s.SetOptFit(11111)
#     s.SetTextColor(ROOT.kBlue)
    htot_DY_0p01.SetMinimum(10e-3)
    htot_DY_0p01.SetMaximum(10e4)
    ps.c.Update()

    fcn_0p05 = ROOT.TF1("fcn_0p05", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn_0p05.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn_0p05.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn_0p05.SetLineColor(ROOT.kRed)
    htot_DY_0p05.Fit(fcn_0p05, 'SREMV')
    htot_DY_0p05.SetMinimum(10e-3)
    htot_DY_0p05.SetMaximum(10e4)
    ps.c.Update()

    fcn_0p1 = ROOT.TF1("fcn_0p1", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn_0p1.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn_0p1.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn_0p1.SetLineColor(ROOT.kGreen)
    htot_DY_0p1.Fit(fcn_0p1, 'SREMV')
    htot_DY_0p1.SetMinimum(10e-3)
    htot_DY_0p1.SetMaximum(10e4)
    ps.c.Update()

    fcn_0p15 = ROOT.TF1("fcn_0p15", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn_0p15.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn_0p15.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn_0p15.SetLineColor(ROOT.kBlack)
    htot_DY_0p15.Fit(fcn_0p15, 'SREMV')
    htot_DY_0p15.SetMinimum(10e-3)
    htot_DY_0p15.SetMaximum(10e4)
    ps.c.Update()

    fcn_0p2 = ROOT.TF1("fcn_0p2", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn_0p2.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn_0p2.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn_0p2.SetLineColor(ROOT.kMagenta)
    htot_DY_0p2.Fit(fcn_0p2, 'SREMV')
    htot_DY_0p2.SetMinimum(10e-3)
    htot_DY_0p2.SetMaximum(10e4)
    ps.c.Update()   

#     s.Draw("e")
#     s_2.Draw("e")

    ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
    
    htot_DY_0p01.GetXaxis().SetRangeUser(600,1200)
    htot_DY_0p01.SetMinimum(1)
    htot_DY_0p01.SetMaximum(500)
    ps.save('mass%i_%i_ZOOM' % (600, 1200), pdf=True, pdf_log=True)
    
    ps.c.Update()
	

#     xax = htot_DY_0p01.GetXaxis()
#     hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-MC)/MC', htot_DY_0p01.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot_DY_0p01.GetNbinsX()+1))
#     for h in [hres]:
# #        h.GetYaxis().SetLabelSize(0.02)
#         h.SetMarkerStyle(2)
#         h.SetStats(0)
#     for i in xrange(1, hres.GetNbinsX()+1):
#         xlo = xax.GetBinLowEdge(i)
#         xhi = xax.GetBinLowEdge(i+1)
#         if xlo >= lo and xhi <= hi:#and xhi < 500:
#             res = fcn.Integral(xlo, xhi)/(xhi-xlo) - htot_DY_0p01.GetBinContent(i)
#             if htot_DY_0p01.GetBinContent(i) > 0:
#                 hres.SetBinContent(i, res/htot_DY_0p01.GetBinContent(i))
#                 hres.SetBinError(i, htot_DY_0p01.GetBinError(i)/htot_DY_0p01.GetBinContent(i))
# #         if xlo >= lo and xhi <= hi and xhi > 500:
# #             res = fcn_2.Integral(xlo, xhi)/(xhi-xlo) - htot_DY_0p01.GetBinContent(i)
# #             if htot_DY_0p01.GetBinContent(i) > 0:
# #                 hres.SetBinContent(i, res/htot_DY_0p01.GetBinContent(i))
# #                 hres.SetBinError(i, htot_DY_0p01.GetBinError(i)/htot_DY_0p01.GetBinContent(i))
#      
#     hres.SetMinimum(-0.25)
#     hres.SetMaximum(0.25)
#     hres.GetXaxis().SetRangeUser(lo, hi)
#     hres.Draw('e')
#     l1 = ROOT.TLine(lo, 0., hi,  0.)
#     l1.Draw()
#     
#     t = ROOT.TPaveLabel(0.15, 0.825, 0.45, 0.925, "K factor: new", 'brNDC')
#     t.SetTextFont(42)
#     t.SetTextSize(0.5)
#     t.SetBorderSize(0)
#     t.SetFillColor(0)
#     t.SetFillStyle(0)
# #     t.Draw()
#     
#     ps.save('res_%i_%i' % (lo, hi), pdf=True, log=False)
#     
#     ps.c.Update()
#        
#        
#     f_data = ROOT.TFile('/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data_OK/NoScale_YesEtaCut_Run2016MuonsOnly/ana_datamc_data.root')
#     c_data = f_data.Our2016MuonsPlusMuonsMinusHistos
#     data = c_data.Get(variable)
#     data.Rebin(rebin)
#     data.SetStats(0)
#     data.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
#     data.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
# #     data.Draw()
#     data.SetLineColor(ROOT.kRed)
#     htot_DY_0p01.SetLineColor(ROOT.kBlue)
#     htot_DY_0p01.GetXaxis().SetRangeUser(lo, hi)
#     htot_DY_0p01.Draw()
#     htot_DY_0p01.Fit(fcn, 'SREMV')
#     data.Draw("same")
#     legend = ROOT.TLegend(0.3, 0.75, 0.6, 0.85)
#     legend.AddEntry(data, "Data", "l")
#     legend.AddEntry(htot_DY_0p01, "MC", "l")
#     legend.AddEntry(fcn, "Fit MC", "l")
#     legend.Draw()
#     htot_DY_0p01.GetXaxis().SetRangeUser(lo, high_for_data)
#     data.GetXaxis().SetRangeUser(lo, high_for_data)
# #     htot_DY_0p01.SetMinimum(10e-5)
# #     htot_DY_0p01.SetMaximum(10e5)
#     data.SetMinimum(10e4)
#     data.SetMaximum(10e5)
#     ps.save('data_%i_%i' % (lo, hi), pdf_log=True)
#     
# 
#     
#     data.SetMinimum(10)
#     data.GetXaxis().SetRangeUser(600,1200)
#     ps.save('data%i_%i_ZOOM' % (600, 1200), pdf_log=True)
#     
#     
#     data.GetXaxis().SetRangeUser(lo, high_for_data)
#     hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-data)/data', data.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(data.GetNbinsX()+1))
#     xax = data.GetXaxis()
#     for h in [hres]:
# #        h.GetYaxis().SetLabelSize(0.02)
#         h.SetMarkerStyle(2)
#         h.SetStats(0)
#     for i in xrange(1, data.GetNbinsX()+1):
#         xlo = xax.GetBinLowEdge(i)
#         xhi = xax.GetBinLowEdge(i+1)
#         if xlo >= lo and xhi <= hi:
# #             print data.GetBinContent(i)
#             res = fcn.Integral(xlo, xhi)/(xhi-xlo) - data.GetBinContent(i)
# #             print res, data.GetBinContent(i), xlo, xhi
#             if data.GetBinContent(i) > 0:
#                 hres.SetBinContent(i, res/data.GetBinContent(i))
#                 hres.SetBinError(i, data.GetBinError(i)/data.GetBinContent(i))
#             if data.GetBinContent(i) == 0:
#                 hres.SetBinContent(i, 1)
#                 hres.SetBinError(i, 0)
# 
# 
#     hres.SetMinimum(-1)
#     hres.SetMaximum( 1)
#     hres.GetXaxis().SetRangeUser(lo, high_for_data)
#     hres.Draw('e')
#     l1 = ROOT.TLine(lo, 0., high_for_data,  0.)
#     l1.Draw()
# 
# #     t = ROOT.TPaveLabel(0.15, 0.825, 0.45, 0.925, "K factor: old", 'brNDC')
# #     t.SetTextFont(42)
# #     t.SetTextSize(0.5)
# #     t.SetBorderSize(0)
# #     t.SetFillColor(0)
# #     t.SetFillStyle(0)
# #     t.Draw()
#     
#     ps.save('res_%i_%i_data' % (lo, hi), pdf=True, log=False)
#  
#     ps.c.Update() 
#  
#  
 
    
#    print fcn.GetProb(), fcn.GetParameter(2), fcn.GetParError(2)
#     c = r.GetCovarianceMatrix()
#     r.Print("V")
#     print " ----------------------------------------------------------------------------------====================== FINE "
#     print c
	
	

#l = range(200, 2200, 200)
# l = [60, 120, 140, 160, 200]
l = [150]
for lo in l:
    print " ----------------------------------------------------------------------------------->>>>>>>>>>>>>>>>>>>>> %d " %lo
    fit_it(lo, high)
    
print " ------------------------------------------------------------------------- Passo al Draw "
print variable, high#, bias
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
# draw_overlay(fsets)
