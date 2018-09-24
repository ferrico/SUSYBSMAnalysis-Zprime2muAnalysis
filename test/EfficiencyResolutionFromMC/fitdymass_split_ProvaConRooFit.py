#!/usr/bin/env python

import os
from ROOT import * 
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)


# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(111)

# variable = 'mass'
variable = 'BB_mass'
# variable = 'BE_mass'

low = fitlow = 150
high = fithigh = 5000
high_for_data = 2500

# ps = plot_saver('plots/fitdymass/split'+ variable)
ps = plot_saver('plots/'+ variable)


int_lumi = 41903.837
rebin = 40
use_non_dy = True

fn = '../DataMCSpectraComparison/MassDistributionForFit.root'
f = ROOT.TFile(fn)
htot = d.Get(variable).Clone()

htot.SetTitle('')
htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.

# bb = ROOT.TFile("./BB_distribution.root", "RECREATE")
# bb.Open()
# htot.Write()
# bb.Close()

htot.Draw()
    
def fit_it(lo, hi):
    
    print " ----------------------------------------------------------------------------------====================== INIZIO "

    htot.Draw()
    ps.c.Update()
    
#     RooRealVar x("x","x",150,5500)
#     RooDataHist MC_mass(htot->GetName(),htot->GetTitle(),RooArgSet(x),RooFit::Import(*htot, kFALSE))

# 	  RooRealVar aL("aL","aL",	24, 19, 29);
# 	  RooRealVar bL("bL","bL",	-2E-3, -2E-2, -2E-4);
# 	  RooRealVar cL("cL","cL",	-1E-7, -7, -2E-7);
# 	  RooRealVar kL("kL","kL",	-3.5, -7, 0);
# 	  RooRealVar aH("aH","aH",	24, 19, 29);
# 	  RooRealVar bH("bH", "bH",	-2E-3, -2E-2, -2E-4);
# 	  RooRealVar cH("cH","cH",	-1E-7, -7, -2E-7);
# 	  RooRealVar dH("dH", "dH",	-1E-10, -1E-11, -1E-9);
# 	  RooRealVar eH("eH", "eH",	-3.5, -7, 0);
# 	  
# 	  RooGenericPdf g("g","g", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
#     					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))",RooArgSet(x,aL,bL,cL,kL, aH, bH, cH, dH, kH)); 
# 	  g.fitTo(MC_mass);
    					   

#     fcn = ROOT.TF1("fcn", "(x<=500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
#     						(x>500)*(exp([4] + [5]*x + [6]*x*x)*x**([7]))", lo, hi)
#     fcn.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "kH")
#     fcn.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -3.5)
    fcn = ROOT.TF1("fcn", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)



#     fcn = ROOT.TF1("fcn", "(x<=500)*(exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])) + \
#     						(x>500)*(exp([5] + [6]*x + [7]*x*x + [8]*x*x*x)*x**([9]))", lo, hi)
#     fcn.SetParNames("aL", "bL", "cL", "dL", "kL", "aH", "bH", "cH", "dH", "kH")
#     fcn.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)

#     fcn = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", lo, hi)
#     fcn = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", 500, hi)
#     fcn = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", lo, 500)
#     fcn.SetParNames("a", "b", "c", "d", "k")
#     fcn.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn.SetLineColor(ROOT.kBlue)
    htot.Fit(fcn, 'SREMV')    
    s = htot.GetListOfFunctions().FindObject("stats")
    s.SetName("Const")
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.69)
    s.SetY2NDC(0.94)
    s.SetOptStat(10)
    s.SetOptFit(11111)
    s.SetTextColor(ROOT.kBlue)


#     fcn_2 = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", 500, hi)
#     fcn_2.SetParNames("a", "b", "c", "d", "k")
#     fcn_2.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5)
#     fcn_2.SetLineColor(ROOT.kRed)
#     htot.Fit(fcn_2, 'SREMV')
#     htot.Fit(fcn, 'SREMV')  
#     s_2 = htot.GetListOfFunctions().FindObject("stats")
#     s_2.SetName("Const")
#     s_2.SetX1NDC(0.73)
#     s_2.SetY1NDC(0.44)
#     s_2.SetY2NDC(0.69)
#     s_2.SetOptStat(10)
#     s_2.SetOptFit(11111)
#     s_2.SetTextColor(ROOT.kRed)


    ps.c.Update()   

#     s.Draw("e")
#     s_2.Draw("e")

    ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
    
    htot.GetXaxis().SetRangeUser(0,500)
    ps.save('mass%i_%i_ZOOM_500' % (lo, hi), pdf=True, pdf_log=True)
    
    htot.GetXaxis().SetRangeUser(0,1000)
    ps.save('mass%i_%i_ZOOM_1000' % (lo, hi), pdf=True, pdf_log=True)
        
    ps.c.Update()
	

    xax = htot.GetXaxis()
    hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-hist)/hist', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
    for h in [hres]:
#        h.GetYaxis().SetLabelSize(0.02)
        h.SetMarkerStyle(2)
        h.SetStats(0)
    for i in xrange(1, hres.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= lo and xhi <= hi:#and xhi < 500:
            res = fcn.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
            if htot.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/htot.GetBinContent(i))
                hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))
#         if xlo >= lo and xhi <= hi and xhi > 500:
#             res = fcn_2.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
#             if htot.GetBinContent(i) > 0:
#                 hres.SetBinContent(i, res/htot.GetBinContent(i))
#                 hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))
     
    hres.SetMinimum(-0.25)
    hres.SetMaximum(0.25)
    hres.GetXaxis().SetRangeUser(lo, hi)
    hres.Draw('e')
    l1 = ROOT.TLine(lo, 0., hi,  0.)
    l1.Draw()
    
    t = ROOT.TPaveLabel(0.15, 0.825, 0.45, 0.925, "K factor: new", 'brNDC')
    t.SetTextFont(42)
    t.SetTextSize(0.5)
    t.SetBorderSize(0)
    t.SetFillColor(0)
    t.SetFillStyle(0)
#     t.Draw()
    
    ps.save('res_%i_%i' % (lo, hi), pdf=True, log=False)
    
    ps.c.Update()
       
#        
#     f_data = ROOT.TFile('ana_datamc_data.root')
#     c_data = f_data.Our2016MuonsPlusMuonsMinusHistos
#     data = c_data.Get(variable)
#     data.Rebin(rebin)
#     data.SetStats(0)
#     data.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
#     data.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
#     data.Draw()
#     htot.SetLineColor(ROOT.kRed)
#     htot.GetXaxis().SetRangeUser(lo, hi)
#     htot.Draw("same")
#     htot.Fit(fcn, 'SREMV')
#     htot.GetXaxis().SetRangeUser(lo, high_for_data)
#     data.GetXaxis().SetRangeUser(lo, high_for_data)
# #     htot.SetMinimum(10e-5)
# #     htot.SetMaximum(10e5)
#     data.SetMinimum(10e-6)
#     data.SetMaximum(10e5)
#     ps.save('data_%i_%i' % (lo, hi), pdf_log=True)
#     
#     data.SetMinimum(10)
#     data.GetXaxis().SetRangeUser(0,500)
#     ps.save('data%i_%i_ZOOM_500' % (lo, hi), pdf_log=True)
# 
#     data.SetMinimum(1)    
#     data.GetXaxis().SetRangeUser(0,1000)
#     ps.save('data%i_%i_ZOOM_1000' % (lo, hi), pdf_log=True)
#     
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

#     t = ROOT.TPaveLabel(0.15, 0.825, 0.45, 0.925, "K factor: old", 'brNDC')
#     t.SetTextFont(42)
#     t.SetTextSize(0.5)
#     t.SetBorderSize(0)
#     t.SetFillColor(0)
#     t.SetFillStyle(0)
#     t.Draw()
    
    ps.save('res_%i_%i_data' % (lo, hi), pdf=True, log=False)
 
    ps.c.Update() 
 
 
 
    
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
print variable, high
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
# draw_overlay(fsets)
