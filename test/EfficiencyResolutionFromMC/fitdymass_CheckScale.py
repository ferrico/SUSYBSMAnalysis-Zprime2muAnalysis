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

# variable = 'DimuonMassVertexConstrained'
# variable = 'DimuonMassVertexConstrained_bb'
# variable = 'DimuonMassVertexConstrained_be'
variable = 'DileptonMass_bb'
# variable = 'DileptonMass_be'
# variable = 'DileptonMass'

# kind = 'MatrixValue'
# kind = 'NoSigmaRestriction'
kind = 'SigmaRestriction'

# sovra = True

if 'Matrix' in kind:
	ps = plot_saver('plots/Matrix/' + variable)
else:
	ps = plot_saver('plots/SigmaRestriction/' + variable)
	
# if '_be' in variable:
# 	category = 'BE + EE'
# 	ps = plot_saver('plots/BE_scale/' + variable)
# elif '_bb' in variable:
# 	category = 'Barrel - Barrel'
# 	ps = plot_saver('plots/BB_resolution/' + variable)
# else:
# 	category = 'Inclusive'
# 	ps = plot_saver('plots/Inclusive_scale/' + variable)
# variable = 'DileptonPt'
# variable = 'LeptonPt'
# variable = 'DileptonMass_bb'



low = fitlow = 150
high = fithigh = 5000
high_for_data = 5000


int_lumi = 36238.734 #36295.39#

# if '_bb' in variable:
# 	rescale_factor = 0.9880
# elif '_be' in variable:
# 	rescale_factor = 0.9625
# else:
# 	rescale_factor = 0.9714

rebin = 80
use_non_dy = True

masses  = ['dy50to120', 'dy120to200', 'dy200to400', 'dy400to800', 'dy800to1400', 'dy1400to2300', 'dy2300to3500', 'dy3500to4500', 'dy4500to6000']
nevents = [2977600, 100000, 100000,  98400, 100000, 95106, 100000, 100000, 100000]
sigmas  = [  1975, 19.32, 2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7]

weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]
# weights = [rescale_factor * int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]

hists_DY = []
hists_dir_DY = '../DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/'
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
if 'Matrix' in kind:
	hists_dir_DY_systematic = '../DataMCSpectraComparison/mc/DY_matrixValue/'
else:
	hists_dir_DY_systematic = '../DataMCSpectraComparison/mc/DY_SigmaRestriction/'

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


hists_dir = '../DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/' 
if use_non_dy:
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
    non_dy_samples = [		WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, 
    						WZ,
#     						WZ_skim,
							ZZ,
# 							ZZ_skim,
    						ZZ_ext,
#     						ZZ_ext_skim,
    						WZ_ext, 
    						Wantitop, tW, 
    						Wjets, 
    						ttbar_lep50to500,
#     						ttbar_lep50to500_OLD, 
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




# if sovra:
# 	htot_DY.Draw()
# 	htot_DY_systematic.Draw("same")
# else:
# 	htot_DY_systematic.Draw()

# htot_DY.SetMinimum(10e2)
# htot_DY_systematic.SetMinimum(10e2)
# htot_DY.SetMaximum(10e4)
# htot_DY.GetXaxis().SetRangeUser(150,500)

legend = ROOT.TLegend(0.25, 0.8, 0.65, 0.95)
legend.AddEntry(htot_DY, "MC", "l")
if 'Matrix' in kind:
	legend.AddEntry(htot_DY_systematic, "MC - DY scale with matrix value", "l")
else:
	legend.AddEntry(htot_DY_systematic, "MC - DY scale #sigma restriction", "l")
legend.Draw()

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
    
    ps.c.Update()
# 
    ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
    
#     htot_DY_0p01.GetXaxis().SetRangeUser(600,1200)
#     htot_DY_0p01.SetMinimum(1)
#     htot_DY_0p01.SetMaximum(500)
#     ps.save('mass%i_%i_ZOOM' % (600, 1200), pdf=True, pdf_log=True)
    
    
    
#     ratio = htot_DY.Clone()
    ratio.Divide(htot_DY_systematic_2)
#     ratio = ratio - 1
    ratio.SetLineColor(1)
    ratio.SetLineWidth(1)
    for i in xrange(1, ratio.GetNbinsX()+1):
    	bin_old = ratio.GetBinContent(i)
    	ratio.SetBinContent(i, bin_old -1)
#     ratio.SetTitle(" ")
    ratio.GetYaxis().SetTitle("Nominal / Scale - 1")
    ratio.SetStats(1)
#     ratio.SetOptStats(111)
    line = ROOT.TLine(70, 0, 5000, 0)
    line.SetLineColor(2)
    line.SetLineWidth(1)
    ratio.Draw()
    ratio.GetXaxis().SetRangeUser(70,hi)
    ps.c.Update()
    print "========================================================================================"    
#     if '_beee' in variable:
    print "======================================================================================="
    fcn_ratio = ROOT.TF1("fcn_ratio", "[0] + [1]*x")
#     fcn_ratio = ROOT.TF1("fcn_ratio", "[0] + [1] * x + [2] * x*x")
# 	    fcn_ratio.SetParNames("m", "q")
# 	    fcn_ratio.SetParameters(1, 4*10e-5) 
    print "ok"
    ratio.Fit(fcn_ratio)#, 'SREMV same')
    sss = ratio.GetListOfFunctions().FindObject("stats")
    sss.SetX1NDC(0.73)
    sss.SetY1NDC(0.69)
    sss.SetY2NDC(0.94)
    sss.SetOptStat(10)
    sss.Draw("same")
    ps.c.Update()
#     	fcn_ratio_2.SetLineColor(ROOT.kBlue)
#     	ratio.Fit(fcn_ratio_2, ' same')
#     	print "ok"
#     	ssss = ratio.GetListOfFunctions().FindObject("stats")
#     	ssss.SetTextColor(ROOT.kBlue)
#     	ssss.SetX1NDC(0.43)
#     	ssss.SetX2NDC(0.73)
#     	ssss.SetY1NDC(0.69)
#     	ssss.SetY2NDC(0.94)
#     	ssss.SetOptStat(10)    	
#     	ssss.Draw("same")
    ratio.GetYaxis().SetRangeUser(-0.5, 0.5)
    line.Draw()

    


#     xax = htot_DY.GetXaxis()
#     hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-MC)/MC', htot_DY.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot_DY.GetNbinsX()+1))
#     for h in [hres]:
# #        h.GetYaxis().SetLabelSize(0.02)
#         h.SetMarkerStyle(2)
#         h.SetStats(0)
#     for i in xrange(1, hres.GetNbinsX()+1):
#         xlo = xax.GetBinLowEdge(i)
#         xhi = xax.GetBinLowEdge(i+1)
#         if xlo >= lo and xhi <= hi:#and xhi < 500:
#             res = fcn.Integral(xlo, xhi)/(xhi-xlo) - fcn_systematic.Integral(xlo, xhi)/(xhi-xlo)
#             if fcn.Integral(xlo, xhi)/(xhi-xlo) > 0:
#                 hres.SetBinContent(i, res/fcn.Integral(xlo, xhi)/(xhi-xlo))
#                 hres.SetBinError(i, htot_DY_0p01.GetBinError(i)/htot_DY_0p01.GetBinContent(i))
#         if xlo >= lo and xhi <= hi and xhi > 500:
#             res = fcn_2.Integral(xlo, xhi)/(xhi-xlo) - htot_DY_0p01.GetBinContent(i)
#             if htot_DY_0p01.GetBinContent(i) > 0:
#                 hres.SetBinContent(i, res/htot_DY_0p01.GetBinContent(i))
#                 hres.SetBinError(i, htot_DY_0p01.GetBinError(i)/htot_DY_0p01.GetBinContent(i))
#      
#     hres.SetMinimum(-0.01)
#     hres.SetMaximum(0.01)
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
    ps.save('res_%i_%i' % (lo, hi), pdf=True, log=False)
#     
    ps.c.Update()
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
print variable, high,int_lumi, sovra, use_non_dy, rebin, kind#, bias
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
# draw_overlay(fsets)
