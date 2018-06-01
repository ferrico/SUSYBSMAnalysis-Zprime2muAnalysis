#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)

# ROOT.gStyle.SetOptStat(0)

variable = 'DimuonMassVertexConstrained'
# variable = 'DimuonMassVertexConstrained_bb'
# variable = 'DimuonMassVertexConstrained_be'

low = fitlow = 120
high = fithigh = 6000

ps = plot_saver('plots/fitdymass/'+ variable)

int_lumi = 37196.77
rebin = 20
use_non_dy = True

# Masses = [50, 120, 200, 400, 800, 1400, 2300, 3500, 4500]


masses  = ['dy50to120', 'dy120to200', 'dy200to400', 'dy400to800', 'dy800to1400', 'dy1400to2300', 'dy2300to3500', 'dy3500to4500', 'dy4500to6000']
nevents = [2977600, 100000, 100000,  98400, 100000, 100000, 100000, 100000, 100000]
#sigmas  = [  1915.,  12.2,  1.53, 0.0462, 0.00586, 0.00194, 1.70e-4, 2.21e-5] # in pb, PYTHIA*1.3
sigmas  = [  1975, 19.32, 2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7]
weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]
#weights = [x/weights[-1] for x in weights]

hists = []
hists_dir = '../DataMCSpectraComparison/mc/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir + fn
    f = ROOT.TFile(fn)
#     d = f.Our2012MuonsPlusMuonsMinusHistos
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    print m, 'w = ', w
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(60, 5500)
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists.append(h)
 
hists_const = []   
hists_dir_const = '../DataMCSpectraComparison/mc/DY_CONSTANT_K/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir_const + fn
    f = ROOT.TFile(fn)
#     d = f.Our2012MuonsPlusMuonsMinusHistos
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    print m, 'w = ', w
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(60, 5500)
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s_Const' % m)
    hists_const.append(h)

if use_non_dy:
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
    non_dy_samples = [WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, 
    						WZ, ZZ,
    						WZ_ext, ZZ_ext,
    						Wantitop, tW, 
    						Wjets, 
#     						ttbar_lep, 
    						ttbar_lep50to500, ttbar_lep_500to800, ttbar_lep_800to1200, ttbar_lep_1200to1800, ttbar_lep1800toInf,
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
#         d = f.Our2012MuonsPlusMuonsMinusHistos
        d = f.Our2016MuonsPlusMuonsMinusHistos

        w = sample.partial_weight * int_lumi
        print sample.name, 'w = ', w
        h = d.Get(variable).Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(60, 5500)
        h.Scale(w)
        h.Draw()
        ps.save('rawmass_%s' % sample.name)
        hists.append(h)
        hists_const.append(h)
  
  
htot = hists[0].Clone()
htot_const = hists_const[0].Clone()

for j in xrange(1, len(hists)):
    htot.Add(hists[j])
    print j
for j in xrange(1, len(hists_const)):
    htot_const.Add(hists_const[j])
    print j

htot.SetTitle('')
htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
htot_const.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
htot_const.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000))
htot.Draw()
htot_const.SetLineColor(ROOT.kRed)
htot_const.Draw("same")

# htot.SetStats(0)
# htot_const.SetStats(0)

#htot.Scale(1/htot.Integral(htot.FindBin(300), htot.FindBin(2000)))

#file = ROOT.TFile('dy.root', 'recreate')
#htot.SetDirectory(file)
#file.Write()
#file.Close
#raise 'done'
    
def fit_it(lo, hi):
#    fcn = ROOT.TF1('fcn', 'exp([0] + [1] * x**[2])', lo, hi)
#    fcn = ROOT.TF1('fcn', 'exp([0] + [1]*x)*x**[2]', lo, hi)
    
    print " ----------------------------------------------------------------------------------====================== INIZIO "

    fcn_const = ROOT.TF1("fcn_const", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", lo, hi)
    fcn_const.SetParNames("a", "b", "c", "d", "k")
    fcn_const.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn_const.SetLineColor(ROOT.kRed)    
    htot_const.Fit(fcn_const, 'SREMV')
#     ROOT.TFitResultPtr 
    
#     s_const = htot_const.GetListOfFunctions().FindObject("stats")
#     s_const.SetTitle("m function")
    s_const = htot_const.FindObject("stats")
    s_const.SetX1NDC(0.63)
    s_const.SetX2NDC(0.93)
    s_const.SetY1NDC(0.74)
    s_const.SetY2NDC(0.99)
    s_const.SetOptStat(0)
    s_const.SetOptFit(11111)
    s_const.SetTextColor(ROOT.kRed)
    

#     ps.c.Update()
    
    fcn = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", lo, hi)
#     fcn = ROOT.TF1("fcn", "exp([0] - [1]*x + [2]*x*x)*x**(-[3])", lo, hi)
    fcn.SetParNames("a", "b", "c", "d", "k")
    fcn.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5)
#     fcn.SetParameters(24, 2E-3, -1E-7, 3.5)
    fcn.SetLineColor(ROOT.kBlue)
    r = htot.Fit(fcn, 'SREMV')    

    s = htot.GetListOfFunctions().FindObject("stats")
#     s.SetTitle("Const")
    s.SetX1NDC(0.33)
    s.SetX2NDC(0.53)
    s.SetY1NDC(0.49)
    s.SetY2NDC(0.74)
    s.SetOptStat(10)
    s.SetOptFit(111)
    s.SetTextColor(ROOT.kBlue)
   

    s_const.Draw("same")
    s.Draw("same")
        
    legend = ROOT.TLegend(0.3,0.85,0.6,0.9)
    legend.AddEntry(htot_const,"k factor: constant", "l")	
    legend.AddEntry(htot, "k factor: m function", "l")
    legend.Draw()



    

    ps.c.Update()
    

    
    ps.save('mass%i_%i' % (lo, hi))



	
    xax = htot.GetXaxis()
    hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-hist)/hist', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
    for h in [hres]:
#        h.GetYaxis().SetLabelSize(0.02)
        h.SetMarkerStyle(2)
        h.SetStats(0)
    for i in xrange(1, hres.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= lo and xhi <= hi:
            res = fcn.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
            if htot.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/htot.GetBinContent(i))
                hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))

    hres.SetMinimum(-1.)
    hres.SetMaximum( 1.)
    hres.GetXaxis().SetRangeUser(lo, hi)
    hres.Draw('e')
    l1 = ROOT.TLine(lo, 0., hi,  0.)
    l1.Draw()

    t = ROOT.TPaveLabel(0.30, 0.825, 0.70, 0.925 , "K factor: m function", 'brNDC')
    t.SetTextFont(42)
    t.SetTextSize(0.5)
    t.SetBorderSize(0)
    t.SetFillColor(0)
    t.SetFillStyle(0)
    t.Draw()

    ps.save('res_%i_%i' % (lo, hi), log=False)
    
    
    
    
    xax = htot_const.GetXaxis()
    hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-hist)/hist', htot_const.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot_const.GetNbinsX()+1))
    for h in [hres]:
#        h.GetYaxis().SetLabelSize(0.02)
        h.SetMarkerStyle(2)
        h.SetStats(0)
    for i in xrange(1, hres.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= lo and xhi <= hi:
            res = fcn_const.Integral(xlo, xhi)/(xhi-xlo) - htot_const.GetBinContent(i)
            if htot_const.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/htot_const.GetBinContent(i))
                hres.SetBinError(i, htot_const.GetBinError(i)/htot_const.GetBinContent(i))

    hres.SetMinimum(-1.)
    hres.SetMaximum( 1.)
    hres.GetXaxis().SetRangeUser(lo, hi)
    hres.Draw('e')
    l1 = ROOT.TLine(lo, 0., hi,  0.)
    l1.Draw()

    t = ROOT.TPaveLabel(0.30, 0.825, 0.70, 0.925, "K factor: constant", 'brNDC')
    t.SetTextFont(42)
    t.SetTextSize(0.5)
    t.SetBorderSize(0)
    t.SetFillColor(0)
    t.SetFillStyle(0)
    t.Draw()
    
    ps.save('res_%i_%i_const' % (lo, hi), log=False)
    
#     print fcn.GetProb(), fcn.GetParameter(2), fcn.GetParError(2)
#     c = r.GetCovarianceMatrix()
#     r.Print("V")
#     print " ----------------------------------------------------------------------------------====================== FINE "
#     print c
	
	
	
def draw_overlay(fsets):
    
    print " ----------------------------------------------------------------------------------====================== INIZIO  DRAW"
    
    ROOT.gStyle.SetOptStat(0);
    ROOT.gStyle.SetOptFit(1111);
    xmin = 120
    xmax = high
    rebin = 20
    htots = []
    for fset in fsets:
        hists = []
        for m,w in zip(masses, weights):
            fn = 'ana_datamc_dy%i%s.root' % (m, fset) if m != 20 else 'ana_datamc_zmumu.root'
            fn = hists_dir + fn
            f = ROOT.TFile(fn)
#             d = f.Our2012MuonsPlusMuonsMinusHistos
            d = f.Our2016MuonsPlusMuonsMinusHistos

            h = d.Get(variable).Clone('dy%i%s' % (m, fset))
            h.Rebin(rebin)
            h.GetXaxis().SetRangeUser(xmin, xmax)
            h.Scale(w)
            hists.append(h)

        htot = hists[0].Clone('htot')
        for j in xrange(1, len(hists)):
            htot.Add(hists[j])
        htot.SetTitle('')
        htots.append(htot)

    for i in xrange(0, len(htots)):
        if i == 0:
            htots[i].Draw()
        else:
            htots[i].Draw("same")
    htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
    htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    ps.save('mass_overlay')

    l = ROOT.TLine(xmin, 1., xmax, 1.)
    for i in xrange(1, len(htots)):
        ps.c.Clear()
        ratio = htots[i].Clone()
        ratio.Divide(htots[0])
        ratio.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
        ratio.GetYaxis().SetTitle('C%i / STARTUP' % i)
        ratio.SetMinimum(0.5)
        ratio.SetMaximum(1.5)
        ratio.SetMarkerStyle(20)
        ratio.Draw("pe")
        ffunc = ROOT.TF1('ffunc', '[0] + x*[1] + x*x*[2]', xmin, xmax)
        ratio.Fit(ffunc, 'IVR')
        l.Draw()
        ps.save('ratio%i_0' % i, log=False)
        if i > 1:
            ratio = htots[i].Clone()
            ratio.Divide(htots[i-1])
            ratio.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
            ratio.GetYaxis().SetTitle('C%i / C%i' % (i, i-1))
            ratio.SetMinimum(0.5)
            ratio.SetMaximum(1.5)
            ratio.SetMarkerStyle(24)
            ratio.SetMarkerColor(ROOT.kBlue)
            ffunc = ROOT.TF1('ffunc', '[0] + x*[1] + x*x*[2]', xmin, xmax)
            ffunc.SetLineColor(ROOT.kBlue)
            ratio.Fit(ffunc, 'IVR')
            ratio.Draw("pe")
            l.Draw()
            ps.save('ratio%i_%i' % (i, i-1), log=False)
        
        print " ----------------------------------------------------------------------------------====================== FINE DRAW "


#l = range(200, 2200, 200)
# l = [60, 120, 140, 160, 200]
l = [120]
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
