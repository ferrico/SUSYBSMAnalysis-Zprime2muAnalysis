#!/usr/bin/env python

##############################
# Histogram making functions #
##############################

def make_mc_hist(do_non_dy,doPI,int_lumi,low,high,rebin,hists_dir,category):
    '''
    - Make total MC histogram
    - inputs are [MCSamples.sample]
    - output is a TH1
    '''
    dy_samples = [dy50to120, dy120to200, dy200to400, dy400to800, dy800to1400, dy1400to2300, dy2300to3500, dy3500to4500, dy4500to6000,dyInclusive50]
    non_dy_samples = [ttbar_lep, tW, Wantitop, WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, WZ, ZZ, Wjets, qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200]#ttbar_pow,
    #non_dy_samples = [ttbar_pow, tW, Wantitop, WWinclusive, WZ, ZZ]
    hists_dy = []
    for sample in dy_samples:
        fn = 'ana_datamc_%s.root' % sample.name
        fn = hists_dir + fn
        f = ROOT.TFile(fn)
        w = sample.partial_weight * int_lumi
        #print sample.name, 'w = ', w
        #h = d.Get('DileptonMass').Clone('%s' % sample.name)
        h = f.Get('Our2012MuonsPlusMuonsMinusHistos').Get(category).Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        hists_dy.append(h)
    hists_nody = []
    for sample in non_dy_samples:
        fn = 'ana_datamc_%s.root' % sample.name
        fn = hists_dir + fn
        f = ROOT.TFile(fn)
        w = sample.partial_weight * int_lumi
        #print sample.name, 'w = ', w
        h = f.Get('Our2012MuonsPlusMuonsMinusHistos').Get(category).Clone('%s' % sample.name)
        #h = d.Get('DileptonMass').Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        hists_nody.append(h)
    hists = []
    # add together DY separately in case of PI background addition
    histDY = hists_dy[0].Clone('histDY')
    for j in xrange(1, len(hists_dy)):
        histDY.Add(hists_dy[j])
    histDY.Draw()
    if doPI:
        histDY = addPI(histDY,rebin,low,high,category,int_lumi)
    # now add all backgrounds
    if do_non_dy:
        hists = hists + hists_nody
        hists.append(histDY)
    else:
        hists.append(histDY)
    histMC = hists[0].Clone('histMC')
    for j in xrange(1, len(hists)):
        histMC.Add(hists[j])
    histMC.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    histMC.GetXaxis().SetLabelFont(42)
    # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    histMC.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin,int_lumi/1000)) 
    histMC.GetYaxis().SetLabelFont(42)
    return histMC

def addPI(histDY,rebin,low,high,category,int_lumi):
    '''
    - Multiply input DY histgram by Photon Induced background
      cross section ratio, R = (DY+PI)/DY
    - Returns TH1
    '''
    PI = ROOT.TF1('pi','[0] + [1]*x + [2]*x*x',low,high)
    #PI.SetParameters(1.044e+00,2.200e-05,1.053e-08) # From Dimitri Method 0 - Paper 2016
    PI.SetParameters(1.025e+00,7.188e-06,1.457e-09) # From Dimitri Method 1 - ICHEP 2016
    histPI = hist_it(PI,rebin,low,high)
    histDYPI = ROOT.TH1F('Hist', ';m(#mu^{+}#mu^{-}) [GeV];Events', 20000, 0, 20000)
    histDYPI.Rebin(rebin)
    histDYPI.GetXaxis().SetRangeUser(low, high)
    histDYPI.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    histDYPI.GetXaxis().SetLabelFont(42)
    # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    histDYPI.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin,int_lumi/1000)) 
    histDYPI.GetYaxis().SetLabelFont(42)
    nBins = histPI.GetNbinsX()
    for i in xrange(1,nBins):
        dy = histDY.GetBinContent(i)
        dye = histDY.GetBinError(i)
        pi = histPI.GetBinContent(i)
        histDYPI.SetBinContent(i,dy*pi)
        histDYPI.SetBinError(i,dye*1.5)
    histDY.SetStats(0)
    histDYPI.SetStats(0)
    histDYPI.Draw()
    histDY.Draw('same')
    histDYPI.SetLineColor(ROOT.kOrange+1)
    histDYPI.SetLineWidth(2)
    histDY.SetLineWidth(2)
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(histDY,'Nominal Background','LE')
    lg.AddEntry(histDYPI,'Nominal+PI Background','LE')
    lg.Draw()
    if category == 'DimuonMassVertexConstrained_bb':
        title='Barrel-Barrel'
        savename = 'BB'
    elif category == 'DimuonMassVertexConstrained_pe':
        title='Barrel-Positive Endcap'
        savename = 'BEp'
    elif category == 'DimuonMassVertexConstrained_bb':
        title='Barrel-Negative Endcap'
        savename = 'BEn'
    else:
        title='All categories'
        savename = 'All'
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    ps.save('%s_PI_hist_overlay'%(savename),pdf=True,log=True)
    return histDYPI

def res_hist(fit,hist,rebin,low,high,resmin=-0.25,resmax=0.25):
    '''
    - Inputs are TF1 fit to a TH1 and the TH1
    - Returns TH1 Residual histogram (fit-hist)/hist
    '''
    hres = ROOT.TH1F('hres', ';m(#mu^{+}#mu^{-}) [GeV];(fit-hist)/hist', 20000, 0, 20000)
    hres.GetXaxis().SetLabelFont(42)
    hres.GetYaxis().SetLabelFont(42)
    hres.Rebin(rebin)
    xax = hres.GetXaxis()
    xax.SetRangeUser(low, high)
    for h in [hres]:
        h.SetMarkerStyle(2)
    for i in xrange(1, hres.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= low and xhi <= high:
            res = fit.Integral(xlo, xhi)/(xhi-xlo) - hist.GetBinContent(i)
            if hist.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/hist.GetBinContent(i))
                hres.SetBinError(i, hist.GetBinError(i)/hist.GetBinContent(i))
    hres.SetMinimum(resmin)
    hres.SetMaximum(resmax)
    return hres

#####################
# Fitting functions #
#####################

def fit_bckg_PI(hist, name,low,high,fitlow, fithigh,rebin,category):
    '''
    - Fit input histogram to background pdf tuned to 
      the DY+PI background shape
    - Returns a TF1
    '''
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])', fitlow, fithigh)
    fit.SetParNames('a', 'b', 'c', 'd', 'k')
    fit.SetParameters(24, -5E-4, -5E-8, -5E-12, -4.5)
    print '\nFit Background Shape with PI',name,'\n'
    hist.Fit(fit,'REM')
    plot_fit(fit,hist,name,low,high,fitlow,fithigh,rebin,category)
    return fit

def fit_bckg(hist, name,low,high,fitlow, fithigh,rebin,category,plot=True):
    '''
    - Fit input histogram to background pdf
    - Returns a TF1
    '''
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])', fitlow, fithigh)
    fit.SetParNames('a', 'b', 'c', 'd', 'k')
    fit.SetParameters(24, -5E-4, -5E-8, -5E-12, -4.5)
    fit.SetParLimits(4,-10,1)
    print '\nFit Background Shape',name,'\n'
    hist.Fit(fit,'REM')
    if plot: plot_fit(fit,hist,name,low,high,fitlow,fithigh,rebin,category)
    return fit

def fit_line(hist,name,fitlow,fithigh):
    '''
    - Fits residual histograms to a straight line
    - Returns a TF1
    '''
    fit = ROOT.TF1('%s'%name, 'pol0', fitlow, fithigh)
    fit.SetParNames('r')
    print '\nFit line to residual plot',name,'\n'
    hist.Fit(fit,'REM')
    fit.SetLineWidth(2)
    fit.SetLineColor(ROOT.kGreen+2)
    return fit

######################
# Plotting functions #
######################

def plot_fit(fit,hist,name,low,high,fitlow,fithigh,rebin,category):
    '''
    - Plotting function
    - Takes in a TF1 and the histogram it was fitted to
      and plots them together. Saves png and root in
      in linear and log formats.
    - Also plots the residual of the nominal fit to the
      total MC histogram.
    '''
    # Plot Fit of Nominal vs. MC histogram
    hist.SetTitle('')
    hist.Draw()
    fit.SetLineColor(ROOT.kBlue)
    ps.c.Update()
    s = hist.GetListOfFunctions().FindObject('stats')
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(10)
    s.SetOptFit(1111)
    s.Draw()
    if category == 'DimuonMassVertexConstrained_bb':
        title='Barrel-Barrel'
    elif category == 'DimuonMassVertexConstrained_pe':
        title='Barrel-Positive Endcap'
    elif category == 'DimuonMassVertexConstrained_ne':
        title='Barrel-Negative Endcap'
    else:
        title='All categories'
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    ps.save('%s'%(name) ,pdf=True,log=True)
    # Plot Residual of Nominal
    hres = res_hist(fit,hist,rebin,low,high)
    line = fit_line(hres,name,fitlow,fithigh)
    plot_res(hres,line,name,category,fitlow,fithigh)

def plot_two(f1,f1Name,hist1,f2,f2Name,hist2,uncName,extra,category,int_lumi):
    '''
    - Plotting function
      takes in 2 TF1s f1,f2 background pdf shapes
    - Returns nothing but saves a pdf, png, and root file in linear
      and log formats
    '''
    ROOT.gStyle.SetOptStat("n");
    # Set Stats for Nominal
    hist1.SetName('%s'%f1Name)
    hist1.Draw()
    ps.c.Update()
    s1 = hist1.GetListOfFunctions().FindObject('stats')
    s1.SetX1NDC(0.73)
    s1.SetY1NDC(0.75)
    s1.SetOptStat(0)
    s1.SetOptFit(1111)
    s1.SetName('%s'%f1Name)
    ps.c.Update()
    hist2.SetName('%s'%f2Name)
    hist2.Draw()
    ps.c.Update()
    s2 = hist2.GetListOfFunctions().FindObject('stats')
    s2.SetX1NDC(0.73)
    s2.SetY1NDC(0.48)
    s2.SetY2NDC(0.73)
    s2.SetOptStat(0)
    s2.SetOptFit(1111)
    s2.SetName('%s'%f2Name)
    ps.c.Update()
    # Draw Minus
    f1.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    f1.GetYaxis().SetLabelFont(42)
    f1.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f1.GetXaxis().SetLabelFont(42)
    f1.SetLineColor(ROOT.kBlack)
    f1.SetTitle('')
    f1.Draw()
    # Draw Nominal
    f2.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f2.SetLineColor(ROOT.kOrange+1)
    f2.SetTitle('')
    s1.Draw('same')
    s2.Draw('same')
    f2.Draw('same')
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(f2,'%s'%f2Name,'L')
    lg.AddEntry(f1,'%s'%f1Name,'L')
    lg.Draw()
    # Draw label
    if category == 'DimuonMassVertexConstrained_bb':
        title='Barrel-Barrel'
        savename = 'BB'
    elif category == 'DimuonMassVertexConstrained_pe':
        title='Barrel-Positive Endcap'
        savename = 'BEp'
    elif category == 'DimuonMassVertexConstrained_ne':
        title='Barrel-Negative Endcap'
        savename = 'BEn'
    else:
        title='All categories'
        savename = 'All'
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Save
    ps.save('%s_%s_%s'%(savename,uncName,extra),pdf=True,log=True)

def plot_three(hist,fitNom,fitPlus,fitMinus,uncName,category,int_lumi):
    '''
    - Plotting function
      takes in 3 TF1s fitNom,fitPlus,fitMinus and plots them
    - Returns nothing but saves a pdf, png, and root file in linear
      and log formats
    '''
    # Settigns fitNom
    fitNom.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    fitNom.GetYaxis().SetLabelFont(42)
    fitNom.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    fitNom.GetXaxis().SetLabelFont(42)
    fitNom.SetLineColor(ROOT.kBlack)
    fitNom.SetTitle('')
    fitNom.GetXaxis().SetRangeUser(200,5500)
    # Settings fitPlus
    fitPlus.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    fitPlus.GetYaxis().SetLabelFont(42)
    fitPlus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    fitPlus.GetXaxis().SetLabelFont(42)
    fitPlus.SetLineColor(ROOT.kBlue)
    fitPlus.SetTitle('')
    fitPlus.GetXaxis().SetRangeUser(200,5500)
    # Settings fitMinus
    fitMinus.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    fitMinus.GetYaxis().SetLabelFont(42)
    fitMinus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    fitMinus.GetXaxis().SetLabelFont(42)
    fitMinus.SetLineColor(ROOT.kOrange+1)
    fitMinus.SetTitle('')
    fitMinus.GetXaxis().SetRangeUser(200,5500)
    # Draw Fits
    fitMinus.Draw('L')
    fitNom.Draw('Lsame')
    fitPlus.Draw('Lsame')
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    uncLabel = uncName+'%'
    lg.AddEntry(fitNom,'Nominal','L')
    lg.AddEntry(fitPlus,'+ %s'%uncLabel,'L')
    lg.AddEntry(fitMinus,'- %s'%uncLabel,'L')
    lg.Draw()
    # Draw label
    if category == 'DimuonMassVertexConstrained_bb':
        title='Barrel-Barrel'
        savename = 'BB'
    elif category == 'DimuonMassVertexConstrained_pe':
        title='Barrel-Positive Endcap'
        savename = 'BEp'
    elif category == 'DimuonMassVertexConstrained_ne':
        title='Barrel-Negative Endcap'
        savename = 'BEn'
    else:
        title='All categories'
        savename = 'All'
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Save
    ps.save('%s_Fit_%spm'%(savename,uncName),pdf=True,log=True)

def plot_res(hres,line,name,category,fitlow,fithigh):
    '''
    - Plot and save residual of (fit-hist)/hist
    '''
    # Draw residual and get stats of fit
    hres.SetTitle('')
    hres.Draw('e')
    line.SetLineColor(ROOT.kGreen+2)
    ps.c.Update()
    s = hres.GetListOfFunctions().FindObject('stats')
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(0)
    s.SetOptFit(1111)
    s.Draw()
    # Draw Label
    if category == 'DimuonMassVertexConstrained_bb':
        title='Barrel-Barrel'
    elif category == 'DimuonMassVertexConstrained_pe':
        title='Barrel-Positive Endcap'
    elif category == 'DimuonMassVertexConstrained_ne':
        title='Barrel-Negative Endcap'
    else:
        title='All categories'
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Draw line at 0
    #l1 = ROOT.TLine(fitlow, 0., fithigh,  0.)
    #l1.SetLineStyle(2)
    #l1.Draw()
    # Save
    ps.save('%s_res'%name, log=False,pdf=True)

def plot_fits(fits,low,high,fitlow,fithigh,rebin,int_lumi,names):
    '''
    - Take in TF1 array of fits; fits = [fit1, fit2, ... , fitN]
      Array of names; names = ['fit1Name','fit2name', ..., 'fitNname']
    - Return a plot of all the fits
    '''
    colors = [ROOT.kOrange+1, ROOT.kBlue, ROOT.kGreen+2, ROOT.kBlack]
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    draw = ''
    for fit,name,color in zip(fits,names,colors):
        if fit.GetName() == 'BB_Fit':
            title='Barrel-Barrel'
        elif fit.GetName() == 'BEp_Fit':
            title='Barrel-Positive Endcap'
        elif fit.GetName() == 'BEn_Fit':
            title='Barrel-Negative Endcap'
        else:
            title='All categories'
        lg.AddEntry(fit,title,'L')
        fit.SetLineColor(color)
        fit.SetTitle('')
        fit.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
        fit.GetXaxis().SetLabelFont(42)
        # assumes original started out with 1 GeV bins, and the xsec is in pb-1. 
        fit.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) 
        fit.Draw(draw)
        draw = 'same'
    lg.Draw()
    ps.save('AllCategoriesFit',pdf=True,log=True)

def plot_norm_fits(fits,low,high,fitlow,fithigh,rebin,int_lumi,names):
    '''
    - Take in TF1 array of fits; fits = [fit1, fit2, ... , fitN]
      Array of names; names = ['fit1Name','fit2name', ..., 'fitNname']
    - Return a plot of all the fits
    '''
    colors = [ROOT.kOrange+1, ROOT.kBlue, ROOT.kGreen+2, ROOT.kBlack]
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    draw = ''
    fitNorms = []
    # Doesn't work unless the list of functions is made before they're drawn?
    for fit in fits:
        # Need to copy the original fit to make everything work
        new = fit.Clone('%s_copy'%fit.GetName())
        fitNorm = make_norm_func(new,fitlow,fithigh)
        fitNorms.append(fitNorm)
    for fit,name,color in zip(fitNorms,names,colors):
        if fit.GetName() == 'BB_Fit_copy_Norm':
            title='Barrel-Barrel'
        elif fit.GetName() == 'BEp_Fit_copy_Norm':
            title='Barrel-Positive Endcap'
        elif fit.GetName() == 'BEn_Fit_copy_Norm':
            title='Barrel-Negative Endcap'
        else:
            title='All categories'
        lg.AddEntry(fit,title,'L')
        fit.SetLineColor(color)
        fit.SetTitle('')
        fit.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
        fit.GetXaxis().SetLabelFont(42)
        fit.GetYaxis().SetTitle('A.U.')
        fit.Draw(draw)
        draw = 'same'
    lg.Draw()
    ps.save('AllCategoriesNormFit',pdf=True,log=True)

def plot_linear_funcs(uncPlus,uncMinus,uncs,colors,low,fithigh):
    '''
    - Takes in and plots linear functions used for bckg shape robustness test
    '''
    # Linear functions
    lg1 = ROOT.TLegend(0.15,0.65,0.45,0.95)
    lg1.SetFillStyle(0)
    lg1.SetTextFont(42)
    lg1.SetBorderSize(0)
    lg2 = ROOT.TLegend(0.15,0.15,0.45,0.45)
    lg2.SetFillStyle(0)
    lg2.SetTextFont(42)
    lg2.SetBorderSize(0)
    draw = ''
    # Reverse the plus so the legend lines up with the lines nicely
    for plus,uncP,colorP,minus,uncM,colorM in zip(reversed(uncPlus),reversed(uncs),reversed(colors),uncMinus,uncs,colors):
        nameP = '+'+uncP+'%'
        nameM = '-'+uncM+'%'
        lg1.AddEntry(plus,nameP,'L')
        lg2.AddEntry(minus,nameM,'L')
        plus.SetTitle('')
        plus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
        plus.GetXaxis().SetRangeUser(100,fithigh)
        plus.GetYaxis().SetTitle('Scaling')
        plus.SetMinimum(0)
        plus.SetMaximum(2)
        plus.SetLineStyle(1)
        plus.SetLineColor(colorP)
        plus.Draw(draw)
        draw = 'same'
        minus.SetTitle('')
        minus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
        minus.GetXaxis().SetRangeUser(100,fithigh)
        minus.GetYaxis().SetTitle('Scaling')
        minus.SetMinimum(0)
        minus.SetMaximum(2)
        minus.SetLineStyle(2)
        minus.SetLineColor(colorM)
        minus.Draw(draw)
    lg1.Draw()
    lg2.Draw()
    ps.save('linear_tests',log=False,pdf=True)

def plot_syst_func(fit, histMC, uncMinus, uncNames,fitlow, fithigh,category,int_lumi,name):
    '''
    - Takes in fit, original MC histogram, lists of +/- shape deformation functions
      and names, and some other inputs
    '''
    # Assume uncs is list of 'minus' functions
    for unc,uncName in zip(uncMinus,uncNames):
        plus, minus = make_pm_func(fit,unc,uncName,fitlow,fithigh,name)
        plot_three(histMC,fit,plus,minus,uncName,category,int_lumi)

#############
# Utilities #
#############

def make_pm_func(fit, unc,uncName,fitlow,fithigh,name):
    '''
    - Takes in nominal background pdf and uncertainty function 
    - Returns TF1s of plus, minus background pdf shapes
    - Need to test to make naming consistent...
      At the moment the pointer name and the TF1 name need to be
      the same. Ex: func = ROOT.TF1('func',...)
    '''
    # Need to copy the original fit to make everything work
    new = fit.Clone('%s_pm_copy'%fit.GetName())
    plus = ROOT.TF1('plus_%s_%s'%(name,uncName),'(2-%s)*%s'%(unc.GetName(),new.GetName()),fitlow,fithigh)
    minus = ROOT.TF1('minus_%s_%s'%(name,uncName),'%s*%s'%(unc.GetName(),new.GetName()),fitlow,fithigh)
    return plus, minus

def make_norm_func(fit,fitlow,fithigh):
    '''
    - Takes in a TF1
    - Returns the same TF1, but normalized to 1
    '''
    NORM = 1/fit.Integral(fitlow,fithigh)
    return ROOT.TF1('%s_Norm'%fit.GetName(),'(%s)*(%s)'%(NORM,fit.GetName()),fitlow,fithigh)

def hist_it(func, rebin,fitlow, fithigh):
    '''
    - Histogram-ize input function
    - Returns a TH1
    '''
    hist = ROOT.TH1F('Hist', ';m(#mu^{+}#mu^{-}) [GeV];Events', 20000, 0, 20000)
    hist.Rebin(rebin)
    hist.GetXaxis().SetRangeUser(fitlow, fithigh)
    for i in xrange(1,hist.GetNbinsX()+1):
        xlo = hist.GetXaxis().GetBinLowEdge(i)
        xhi = hist.GetXaxis().GetBinLowEdge(i+1)
        if xlo >= fitlow and xhi <= fithigh:
            integ = func.Integral(xlo,xhi) / (xhi - xlo) 
            hist.SetBinContent(i, integ)
    return hist

def makeKfactorPlot(kFactor,nominal,fitlow,fithigh,int_lumi,name):
    '''
    - asdf
    - asdf
    '''
    # Copy nominal fit and set styles
    nomCopy = nominal.Clone('nomCopy_%s'%name)
    nomCopy.SetTitle('')
    nomCopy.SetLineColor(ROOT.kBlue)
    # New background shape = kfactor function * nominal background shape
    # Set styles
    kFactorBckgShape = ROOT.TF1('kFactorBckgShape','nomCopy_%s*kFactor'%(name),fitlow,fithigh)
    kFactorBckgShape.SetTitle('')
    kFactorBckgShape.SetLineColor(ROOT.kOrange+1)
    kFactorBckgShape.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) 
    kFactorBckgShape.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    # Draw Functions
    kFactorBckgShape.Draw()
    nomCopy.Draw('same')
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(nomCopy,'Nominal','L')
    lg.AddEntry(kFactorBckgShape,'k-factor applied','L')
    lg.Draw()
    # Draw title
    if name=='BB': title = 'Barrel-Barrel'
    elif name=='BEp': title = 'Barrel - Positive Endcap'
    elif name=='BEn': title = 'Barrel - Negative Endcap'
    else: title = 'All Categories'
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Save
    ps.save('%s_kFactor'%name,log=True)

def plotKfactor(kFactor,fitlow,fithigh):
    '''
    - asdf
    - asdf
    '''
    # kfactor function
    # Set styles
    kFactor.SetTitle('')
    kFactor.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    # Draw Function
    kFactor.Draw()
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(kFactor,'k-factor','L')
    lg.Draw()
    # Save
    ps.save('kFactor',log=True)

#**************************#
# Main portion of the code #
#**************************#

if __name__=='__main__':
    import numpy as np
    import ROOT,sys
    ROOT.gROOT.SetBatch(True)
    from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
    import argparse
    parser = argparse.ArgumentParser(description='submits limits for a given mass')
    parser.add_argument('--tag',help='Tagged directory name for output plots',default='')
    parser.add_argument('--doPI',help='Make Photon-Induced robustness test plots',default=False)
    parser.add_argument('--doShapeSyst',help='Make background shape robustness test plots',default=False)
    parser.add_argument('--doKfactor',help='Factor in k-factor on the background shape',default=False)
    args = parser.parse_args()

    set_zp2mu_style()
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadRightMargin(0.02)
    ps = plot_saver('plots/'+args.tag)
    ROOT.gPad.SetTicks(1)
    ROOT.TH1.AddDirectory(0)

    #hists_dir = '/afs/cern.ch/work/c/cschnaib/Zprime2muAnalysis/CMSSW_8_0_3/DataMCSpectraComparison/mc/80X_v1/'
    # Trigger weight applied
    #hists_dir = '/afs/cern.ch/work/f/ferrico/public/Root_ZPrime_DONT_DELETE/DataMCSpectra/mc/'
    # FIXME set hists_dir to where MC files are
    # Trigger weight down applied
    hists_dir = '/afs/cern.ch/work/f/ferrico/private/Codice_ZPrime_8_NoTrigger/CMSSW_8_0_3_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/mc_TriggerScale/'
    rebin = 40
    low = 40
    fitlow = 200
    fithigh = 5500
    high = 5500
    int_lumi = 1000.

    fits = []
    names = []

    if args.doShapeSyst:
        #unc05 = ROOT.TF1('unc05','[0] + [1]*x',low,fithigh)                      
        #unc05.SetParameters(1.0012,-1.71E-5)
        #unc10 = ROOT.TF1('unc10','[0] + [1]*x',low,fithigh)                      
        #unc10.SetParameters(1.0024,-3.41E-5)
        unc15minus = ROOT.TF1('unc15minus','[0] + [1]*x',low,fithigh)
        unc15minus.SetParameters(1.0036,-5.12E-5)
        unc20minus = ROOT.TF1('unc20minus','[0] + [1]*x',low,fithigh)
        unc20minus.SetParameters(1.0048,-6.83E-5)
        unc25minus = ROOT.TF1('unc25minus','[0] + [1]*x',low,fithigh)
        unc25minus.SetParameters(1.0060,-8.53E-5)
        #unc30minus = ROOT.TF1('unc30','[0] + [1]*x',low,fithigh)                      
        #unc30minus.SetParameters(1.0072,-10.24E-5)
        #unc35minus = ROOT.TF1('unc35','[0] + [1]*x',low,fithigh)                      
        #unc35minus.SetParameters(1.0084,-11.94E-5)
        #unc40minus = ROOT.TF1('unc40','[0] + [1]*x',low,fithigh)                      
        #unc40minus.SetParameters(1.0096,-13.65E-5)
        #unc45minus = ROOT.TF1('unc40','[0] + [1]*x',low,fithigh)                      
        #unc45minus.SetParameters(1.0107,-15.36E-5)
        unc50minus = ROOT.TF1('unc50minus','[0] + [1]*x',low,fithigh)
        unc50minus.SetParameters(1.0119,-17.06E-5)
        unc15plus = ROOT.TF1('unc15plus','2-unc15minus',low,fithigh)
        unc20plus = ROOT.TF1('unc20plus','2-unc20minus',low,fithigh)
        unc25plus = ROOT.TF1('unc25plus','2-unc25minus',low,fithigh)
        unc50plus = ROOT.TF1('unc50plus','2-unc50minus',low,fithigh)
        # Lists
        uncPlus =  [ unc15plus, unc20plus, unc25plus, unc50plus]
        uncMinus = [unc15minus,unc20minus,unc25minus,unc50minus]
        uncNames = ['15','20','25','50']
        colors =  [ROOT.kBlack,ROOT.kBlue,ROOT.kGreen+2,ROOT.kOrange+1]
        plot_linear_funcs(uncPlus,uncMinus,uncNames,colors,low,fithigh)

    for category in ['DimuonMassVertexConstrained_bb',
                     'DimuonMassVertexConstrained_ne',
                     'DimuonMassVertexConstrained_pe',
                     'DimuonMassVertexConstrained']:
        print category
        if category=='DimuonMassVertexConstrained_bb':
            name = 'BB'
        if category=='DimuonMassVertexConstrained_pe':
            name = 'BEp'
        if category=='DimuonMassVertexConstrained_ne':
            name = 'BEn'
        if category=='DimuonMassVertexConstrained':
            name = 'All'

        # Make nominal MC histogram
        histMC = make_mc_hist(True,False,int_lumi,low,high,rebin,hists_dir,category)
        # Make nominal MC fit and saves plot
        nominal = fit_bckg(histMC,name+'_Fit',low,high,fitlow,fithigh,rebin,category)
        fits.append(nominal)
        names.append(category)

        # Photon-Induced block
        if args.doPI:
            # Make MC histogram with PI background
            histMCPI = make_mc_hist(True,True,int_lumi,low,high,rebin,hists_dir,category)
            nominalPI = fit_bckg_PI(histMCPI,name+'_PI_Fit',low,high,fitlow,fithigh,rebin,category)
            # Plot Bckg and Bckg+PI fit functions
            plot_two(nominal,'Nominal',histMC,nominalPI,'Nominal+PI',histMCPI,'nominal_PI','Fit',category,int_lumi)

        # Plot +/nominal/- background shapes
        if args.doShapeSyst: plot_syst_func(nominal,histMC,uncMinus,uncNames,fitlow,fithigh,category,int_lumi,name)
        
        # Factor in k-factor parametrization from Sam
        if args.doKfactor: 
            #kFactorFunc = ROOT.TF1('kFactorFunc','1.01696-7.73522E-5*x+6.69239E-9*x*x',
            #std::min(1.01696-7.73522E-5*mass+6.69239E-9*mass*mass,1)
            kFactor = ROOT.TF1('kFactor','TMath::Min(1.01696 - 7.73522E-5*x + 6.69239E-9*x*x,1.0)',low,high)
            makeKfactorPlot(kFactor,nominal,fitlow,fithigh,int_lumi,name)
            pass

        print '\n*****************\n'
    print '\n'

    # Plot All categories on top of each other
    plot_fits(fits,low,high,fitlow,fithigh,rebin,int_lumi,names)
    plot_norm_fits(fits,low,high,fitlow,fithigh,rebin,int_lumi,names)
    # Plot k-factor function
    if args.doKfactor: plotKfactor(kFactor,fitlow,fithigh)
    print 'done'


