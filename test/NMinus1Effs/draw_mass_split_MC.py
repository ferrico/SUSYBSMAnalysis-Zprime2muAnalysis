#!/usr/bin/env python

# (py draw.py >! plotay Filippo/nminus1effs/out.draw) && tlp plots/nminus1effs

# Import nm1entry.py
# - from MCSamples and roottools import *
# - nm1entry class
# - nminus1s list
# - pretty dictionary
# - styles dictionary
from nm1entry import *

def draw_mass_test(tag, printStats, lumi, do_tight=False):

    # Specific stylings for this script on top of zp2mu_style
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadRightMargin(0.02)
    ROOT.gStyle.SetTitleX(0.12)
    #ROOT.gStyle.SetTitleH(0.07)
    ROOT.TH1.AddDirectory(0)
    outfile = ROOT.TFile("whargl.root","recreate")
    iarp=0

    # 'plots' = '/afs/cern.ch/work/c/cschnaib/Zprime2muAnalysis/NMinus1Effs/plots/TAG/tag'
    # TAG = MC/Data production TAG
    # tag = sub tag for each version of plots made with a single production TAG
#     psn = 'plots/' + categoria_1 + cartella + '/mass_split/'
    psn = 'plots/'
    if tag:
        psn = 'plots/'

    # cjsbad - fix the referencing to tightnm1...
    #if do_tight:
    #    psn += '_tight'
    #    nminus1s = tightnm1[:]

    ps = plot_saver(psn, size=(600,600), log=False, pdf=True)
    ps.c.SetBottomMargin(0.2)

    data = nm1entry('data', True, lumi)#lumiCD )

    refXS = dy50to120.cross_section
    refN = dy50to120.nevents

    samples = [ dyInclusive50,
# 				qcd80to120, qcd120to170, qcd170to300, qcd300to470, qcd470to600, qcd600to800, 
# 				qcd800to1000, qcd1000to1400, qcd1400to1800, qcd1800to2400, qcd2400to3200, qcd3200,
				Wantitop, tW, 
				ttbar_lep50to500, ttbar_lep_500to800, ttbar_lep_800to1200, ttbar_lep_1200to1800, ttbar_lep_1800,
				Wjets, 
				WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, 
				ZZ, WZ, 
				ZZ_ext, WZ_ext, 
				dy50to120, 
				dy120to200, dy200to400, dy400to800, dy800to1400, 
				dy1400to2300, dy2300to3500, dy3500to4500, dy4500to6000
				]

    samples_DY = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000, dyInclusive50]

    samples_diboson = [WZ,WZ_ext, ZZ, ZZ_ext, WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500]

    samples_Wjets = [ Wjets, 
    				Wantitop,tW, 
    				# qcd80to120, qcd120to170, qcd170to300, qcd300to470, qcd470to600,qcd600to800,
# qcd800to1000, qcd1000to1400, qcd1400to1800, qcd1800to2400, qcd2400to3200,qcd3200
]

    samples_tt = [ttbar_lep50to500, ttbar_lep_500to800, ttbar_lep_800to1200, ttbar_lep_1200to1800, ttbar_lep_1800]

    # All MC samples
    # lumi
    mc_samples = [nm1entry(sample,False,lumi) for sample in samples]
    #for mc_sample in mc_samples:
    #    exec '%s = mc_sample' % mc_sample.name

    mc_samples_DY = [nm1entry(sample,False,lumi) for sample in samples_DY]
#    for mc_sample_DY in mc_samples_DY:
#        exec '%s = mc_sample_DY' % mc_sample_DY.name

    mc_samples_tt = [nm1entry(sample,False,lumi) for sample in samples_tt]
#    for mc_sample_tt in mc_samples_tt:
#        exec '%s = mc_sample_tt' % mc_sample_tt.name

    mc_samples_diboson = [nm1entry(sample,False,lumi) for sample in samples_diboson]
#    for mc_sample_diboson in mc_samples_diboson:
#        exec '%s = mc_sample_diboson' % mc_sample_diboson.name

    mc_samples_Wjets = [nm1entry(sample,False,lumi) for sample in samples_Wjets]
#    for mc_sample_Wjets in mc_samples_Wjets:
#        exec '%s = mc_sample_Wjets' % mc_sample_Wjets.name



    #bin_width = 20
    #maxX = 2500
    #minX = 60
    #nBins = (maxX-minX)/bin_width
    #mass_range = []
    #for i in range(3,nBins):
    #    ibin = i*bin_width
    #    mass_range.append(ibin)
    #print mass_range
#    mass_range = [60,120,180,240,320,500,1000,2500]
#     mass_range = [50,120,200,300,400,500,600,800,1400,2300]
    mass_range = [50,120,200,400,800,1400,2300]#,3500]

    to_use = {
    #   'sample':[MC,Data],
        'NoPt':[mc_samples,data],
        'NoDB':[mc_samples,data],
        'NoIso':[mc_samples,data],
        'NoTkLayers':[mc_samples,data],
        'NoPxHits':[mc_samples,data],
        'NoMuHits':[mc_samples,data],
        'NoMuMatch':[mc_samples,data],
        'NoVtxProb':[mc_samples,data],
        'NoB2B':[mc_samples,data],
        'NoDptPt':[mc_samples,data],
        'NoTrgMtch':[mc_samples,data],
        }
 
    to_use_DY = {
        'NoPt':[mc_samples_DY, data],
        'NoDB':[mc_samples_DY,data],
        'NoIso':[mc_samples_DY,data],
        'NoTkLayers':[mc_samples_DY,data],
        'NoPxHits':[mc_samples_DY,data],
        'NoMuHits':[mc_samples_DY,data],
        'NoMuMatch':[mc_samples_DY,data],
        'NoVtxProb':[mc_samples_DY,data],
        'NoB2B':[mc_samples_DY,data],
        'NoDptPt':[mc_samples_DY,data],
        'NoTrgMtch':[mc_samples_DY,data],
                }
    
    to_use_tt = {
    'NoPt':[mc_samples_tt, data],
    'NoDB':[mc_samples_tt,data],
    'NoIso':[mc_samples_tt,data],
    'NoTkLayers':[mc_samples_tt,data],
    'NoPxHits':[mc_samples_tt,data],
    'NoMuHits':[mc_samples_tt,data],
    'NoMuMatch':[mc_samples_tt,data],
    'NoVtxProb':[mc_samples_tt,data],
    'NoB2B':[mc_samples_tt,data],
    'NoDptPt':[mc_samples_tt,data],
    'NoTrgMtch':[mc_samples_tt,data],
    }
    
    to_use_diboson = {
    'NoPt':[mc_samples_diboson, data],
    'NoDB':[mc_samples_diboson,data],
    'NoIso':[mc_samples_diboson,data],
    'NoTkLayers':[mc_samples_diboson,data],
    'NoPxHits':[mc_samples_diboson,data],
    'NoMuHits':[mc_samples_diboson,data],
    'NoMuMatch':[mc_samples_diboson,data],
    'NoVtxProb':[mc_samples_diboson,data],
    'NoB2B':[mc_samples_diboson,data],
    'NoDptPt':[mc_samples_diboson,data],
    'NoTrgMtch':[mc_samples_diboson,data],
    }
    
    to_use_Wjets = {
    'NoPt':[mc_samples_Wjets, data],
    'NoDB':[mc_samples_Wjets,data],
    'NoIso':[mc_samples_Wjets,data],
    'NoTkLayers':[mc_samples_Wjets,data],
    'NoPxHits':[mc_samples_Wjets,data],
    'NoMuHits':[mc_samples_Wjets,data],
    'NoMuMatch':[mc_samples_Wjets,data],
    'NoVtxProb':[mc_samples_Wjets,data],
    'NoB2B':[mc_samples_Wjets,data],
    'NoDptPt':[mc_samples_Wjets,data],
    'NoTrgMtch':[mc_samples_Wjets,data],
    }     
     

    yrange = {
    #   'sample':    (ymin,ymax),
        'NoPt':      (0.00,1.01),
        'NoDB':      (0.95,1.001),
        'NoIso':     (0.8,1.01),
        'NoTkLayers':(0.8,1.01),
        'NoPxHits':  (0.80,1.01),
        'NoMuHits':  (0.80,1.01),
        'NoMuMatch': (0.1,1.01),
        'NoVtxProb': (0.9,1.01),
        'NoB2B':     (0.80,1.001),
        'NoDptPt':   (0.8,1.01),
        'NoTrgMtch': (0.90,1.001),
        }
    #global_ymin = 0.
    global_ymin = None


    ROOT.gStyle.SetTitleX(0.25)
    ROOT.gStyle.SetTitleY(0.50)

    for nminus1 in nminus1s:
        pretty_name = pretty[nminus1]
        lg = ROOT.TLegend(0.25, 0.21, 0.91, 0.44)
        lg.SetTextSize(0.03)
        lg.SetFillColor(0)
        lg.SetBorderSize(1)
        same = 'A'
        effs = []


##### all simulation
        print 'All simulation'
        for entry in to_use[nminus1]: #,mass_range 
            l = len(mass_range)-1
            nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
            nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))
            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                if printStats: table(entry,nminus1,mass_range)
                color, fill = styles[entry.name]
                hnum = entry.histos['NoNo']
                hden = entry.histos[nminus1]
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                    nminus1_num.SetBinContent(mbin+1, num)
                    nminus1_den.SetBinContent(mbin+1, den)
                eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
            else:
                if printStats: table_wald(mc,nminus1,mass_range)
                p_hats = []
                errsW = []
                x = []
                ex = []
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    numTot = 0
                    denTot = 0
                    err2sum = 0
                    numsW = []
                    densW = []
                    err2s = []
                    for i,mc in enumerate(entry):
                        hnum = mc.histos['NoNo']
                        hden = mc.histos[nminus1]
                        mlow = mass_range[mbin]
                        mhigh = mass_range[mbin+1]
                        numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                        denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                        if numInt<=denInt and denInt>0:
                            p_hat_mc = float(numInt)/denInt
                            err2 = p_hat_mc*(1-p_hat_mc)/denInt
                        elif numInt>denInt:
                            print "numInt>denInt"
                            print "NM1, sample, mlow, mhigh, num, den"
                            print nminus1, mc.name, mlow, mhigh, numInt, denInt
                            print
                            p_hat_mc = 1
                            err2 = 0
                        else:
                            if numInt!=0 or denInt!=0:
                                print "ELSE"
                                print "NM1, sample, mlow, mhigh, num, den"
                                print nminus1, mc.name, mlow, mhigh, numInt, denInt
                                print
                            p_hat_mc = 0
                            err2 = 0
                        numTot = numTot + numInt*mc.partial_weight
                        denTot = denTot + denInt*mc.partial_weight
                        numsW.append(numInt*mc.partial_weight)
                        densW.append(denInt*mc.partial_weight)
                        err2s.append(err2)
                    p_hat = float(numTot)/denTot
                    err2sum = sum( ((m/denTot)**2 * e2) for m,e2 in zip(densW,err2s))
                    if err2sum<0:
                        err2sum = 0
                    elif p_hat==0:
                        err2sum = 0
                    p_hats.append(p_hat)
                    errsW.append(err2sum**0.5)
                    x.append(nminus1_num.GetXaxis().GetBinCenter(mbin+1))
                    ex.append(nminus1_num.GetXaxis().GetBinWidth(mbin+1)/2)
                eff = ROOT.TGraphAsymmErrors(len(x), *[array('d',obj) for obj in (x,p_hats,ex,ex,errsW,errsW)])
            eff.SetTitle(pretty_name)
            ymin, ymax = yrange[nminus1]
            eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin, ymax)
            eff.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
            eff.GetYaxis().SetLabelSize(0.027)
            eff.GetYaxis().SetTitle('n-1 efficiency')
            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                draw = 'P'
                eff.SetLineColor(color)
                eff.SetMarkerStyle(20)
                eff.SetMarkerSize(1.05)
                eff.SetMarkerColor(color)
#                 lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (lumi/1000.), 'LP')
                lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (entry.lumi/1000.), 'LP')
            else:
                if len(entry)==len(mc_samples):
                    draw = '2'
                    eff.SetLineColor(ROOT.kGreen+2)
                    eff.SetFillColor(ROOT.kGreen+2)
                    eff.SetFillStyle(1001)
                    lg.AddEntry(eff,'All simulation','LF')
                else:
                    raise ValueError("Set line+fill color, fill Style, and set legend entry for MC entry list")
            #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
            draw += same
            eff.Draw(draw)
            effs.append(eff)
            same = ' same'
            outfile.cd()
            eff.Write("arp%d"%iarp)
            iarp+=1
        # end for entry in to_use[name]: # entry is a specific sample

###### DY
        #DY
        print 'DY'
        for entry in to_use_DY[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)
       
            l = len(mass_range)-1
            nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
            nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))

            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
            	    if 'data' in entry.name:
            		    continue	
        	#table_wald(entry,nminus1,mass_range)
            	    color, fill = styles[entry.name]
            	    hnum = entry.histos['NoNo']
            	    hden = entry.histos[nminus1]
        	    for mbin in range(len(mass_range)):
        		    if mbin == (len(mass_range)-1): continue
        		    mlow = mass_range[mbin]
        		    mhigh = mass_range[mbin+1]
        		    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
        		    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
        		    nminus1_num.SetBinContent(mbin+1, num)
        		    nminus1_den.SetBinContent(mbin+1, den)
        	    eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
            else:
#             for a,mc in enumerate(entry):
#                 table_wald(mc,nminus1,mass_range)
                p_hats_DY = []
                errsW_DY = []
                x_DY = []
                ex_DY = []
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    numTot = 0
                    denTot = 0
                    err2sum = 0
                    numsW_DY = []
                    densW_DY = []
                    err2s_DY = []
                    for i,mc_DY in enumerate(entry):
                        hnum = mc_DY.histos['NoNo']
                        hden = mc_DY.histos[nminus1]
                        mlow = mass_range[mbin]
                        mhigh = mass_range[mbin+1]
                        numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                        denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                        if numInt<=denInt and denInt>0:
                            p_hat_mc = float(numInt)/denInt
                            err2 = p_hat_mc*(1-p_hat_mc)/denInt
                        elif numInt>denInt:
                            print "numInt>denInt"
                            print "NM1, sample, mlow, mhigh, num, den"
                            print nminus1, mc_DY.name, mlow, mhigh, numInt, denInt
                            print
                            p_hat_mc = 1
                            err2 = 0
                        else:
                            if numInt!=0 or denInt!=0:
                                print "ELSE"
                                print "NM1, sample, mlow, mhigh, num, den"
                                print nminus1, mc_DY.name, mlow, mhigh, numInt, denInt
                                print
                            p_hat_mc = 0
                            err2 = 0
                        numTot = numTot + numInt*mc_DY.partial_weight
#                    print '%s ----- %s ---- %s' % (denTot, denInt, mc.partial_weight)
                        denTot = denTot + denInt*mc_DY.partial_weight
                        numsW_DY.append(numInt*mc_DY.partial_weight)
                        densW_DY.append(denInt*mc_DY.partial_weight)
                        err2s_DY.append(err2)
                    p_hat = float(numTot)/denTot
                    err2sum = sum( ((m/denTot)**2 * e2) for m,e2 in zip(densW_DY,err2s_DY))
                    if err2sum<0:
                        err2sum = 0
                    elif p_hat==0:
                        err2sum = 0
                    p_hats_DY.append(p_hat)
                    errsW_DY.append(err2sum**0.5)
                    x_DY.append(nminus1_num.GetXaxis().GetBinCenter(mbin+1))
                    ex_DY.append(nminus1_num.GetXaxis().GetBinWidth(mbin+1)/2)
                eff = ROOT.TGraphAsymmErrors(len(x_DY), *[array('d',obj) for obj in (x_DY,p_hats_DY,ex_DY,ex_DY,errsW_DY,errsW_DY)])
            eff.SetTitle(pretty_name)
            ymin, ymax = yrange[nminus1]
            eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin, ymax)
            eff.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
            eff.GetYaxis().SetLabelSize(0.027)
            eff.GetYaxis().SetTitle('n-1 efficiency')
            draw = '2'
            eff.SetLineColor(ROOT.kMagenta+2)
            eff.SetFillColor(ROOT.kMagenta+2)
            eff.SetFillStyle(3001)
            lg.AddEntry(eff,'DY','LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
            draw += same
            eff.Draw(draw)
            effs.append(eff)
            same = ' same'
            outfile.cd()
            eff.Write("arp%d"%iarp)
            iarp+=1

##### Diboson
        print 'Diboson'
        for entry in to_use_diboson[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)
       
            l = len(mass_range)-1
            nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
            nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))

            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
        	    if 'data' in entry.name:
        		    continue	
        	#table_wald(entry,nminus1,mass_range)
        	    color, fill = styles[entry.name]
        	    hnum = entry.histos['NoNo']
        	    hden = entry.histos[nminus1]
        	    for mbin in range(len(mass_range)):
        		    if mbin == (len(mass_range)-1): continue
        		    mlow = mass_range[mbin]
        		    mhigh = mass_range[mbin+1]
        		    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
        		    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
        		    nminus1_num.SetBinContent(mbin+1, num)
        		    nminus1_den.SetBinContent(mbin+1, den)
        	    eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
            else:
#             for a,mc in enumerate(entry):
#                 table_wald(mc,nminus1,mass_range)
                p_hats_diboson = []
                errsW_diboson = []
                x_diboson = []
                ex_diboson = []
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    numTot = 0
                    denTot = 0
                    err2sum = 0
                    numsW_diboson = []
                    densW_diboson = []
                    err2s_diboson = []
                    for i,mc_diboson in enumerate(entry):
                        hnum = mc_diboson.histos['NoNo']
                        hden = mc_diboson.histos[nminus1]
                        mlow = mass_range[mbin]
                        mhigh = mass_range[mbin+1]
                        numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                        denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                        if numInt<=denInt and denInt>0:
                            p_hat_mc = float(numInt)/denInt
                            err2 = p_hat_mc*(1-p_hat_mc)/denInt
                        elif numInt>denInt:
                            print "numInt>denInt"
                            print "NM1, sample, mlow, mhigh, num, den"
                            print nminus1, mc_diboson.name, mlow, mhigh, numInt, denInt
                            print
                            p_hat_mc = 1
                            err2 = 0
                        else:
                            if numInt!=0 or denInt!=0:
                                print "ELSE"
                                print "NM1, sample, mlow, mhigh, num, den"
                                print nminus1, mc_diboson.name, mlow, mhigh, numInt, denInt
                                print
                            p_hat_mc = 0
                            err2 = 0
                        numTot = numTot + numInt*mc_diboson.partial_weight
#                    print '%s ----- %s ---- %s' % (denTot, denInt, mc.partial_weight)
                        denTot = denTot + denInt*mc_diboson.partial_weight
                        numsW_diboson.append(numInt*mc_diboson.partial_weight)
                        densW_diboson.append(denInt*mc_diboson.partial_weight)
                        err2s_diboson.append(err2)                
                    if denTot == 0:
                	    p_hat = 0
                	    err2sum = 0
                    else:
	                    p_hat = float(numTot)/denTot
                	    err2sum = sum( ((m/denTot)**2 * e2) for m,e2 in zip(densW_diboson,err2s_diboson))
                    if err2sum<0:
                        err2sum = 0
                    elif p_hat==0:
                        err2sum = 0
                    p_hats_diboson.append(p_hat)
                    errsW_diboson.append(err2sum**0.5)
                    x_diboson.append(nminus1_num.GetXaxis().GetBinCenter(mbin+1))
                    ex_diboson.append(nminus1_num.GetXaxis().GetBinWidth(mbin+1)/2)
                eff = ROOT.TGraphAsymmErrors(len(x_diboson), *[array('d',obj) for obj in (x_diboson,p_hats_diboson,ex_diboson,ex_diboson,errsW_diboson,errsW_diboson)])
            eff.SetTitle(pretty_name)
            ymin, ymax = yrange[nminus1]
            eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin, ymax)
            eff.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
            eff.GetYaxis().SetLabelSize(0.027)
            eff.GetYaxis().SetTitle('n-1 efficiency')
            draw = '2'
            eff.SetLineColor(ROOT.kBlue+2)
            eff.SetFillColor(ROOT.kBlue+2)
            eff.SetFillStyle(3001)
            lg.AddEntry(eff,'Diboson','LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
            draw += same
            eff.Draw(draw)
            effs.append(eff)
            same = ' same'
            outfile.cd()
            eff.Write("arp%d"%iarp)
            iarp+=1                         



#### W jets
        print 'Wjets'
        for entry in to_use_Wjets[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)
       
            l = len(mass_range)-1
            nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
            nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))

            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
        	    if 'data' in entry.name:
        		    continue	
        	#table_wald(entry,nminus1,mass_range)
        	    color, fill = styles[entry.name]
        	    hnum = entry.histos['NoNo']
        	    hden = entry.histos[nminus1]
        	    for mbin in range(len(mass_range)):
        		    if mbin == (len(mass_range)-1): continue
        		    mlow = mass_range[mbin]
        		    mhigh = mass_range[mbin+1]
        		    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
        		    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
        		    nminus1_num.SetBinContent(mbin+1, num)
        		    nminus1_den.SetBinContent(mbin+1, den)
        	    eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
            else:
#             for a,mc in enumerate(entry):
#                 table_wald(mc,nminus1,mass_range)
                p_hats_Wjets = []
                errsW_Wjets = []
                x_Wjets = []
                ex_Wjets = []
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    numTot = 0
                    denTot = 0
                    err2sum = 0
                    numsW_Wjets = []
                    densW_Wjets = []
                    err2s_Wjets = []
                    for i,mc_Wjets in enumerate(entry):
                        hnum = mc_Wjets.histos['NoNo']
                        hden = mc_Wjets.histos[nminus1]
                        mlow = mass_range[mbin]
                        mhigh = mass_range[mbin+1]
                        numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                        denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                        if numInt<=denInt and denInt>0:
                            p_hat_mc = float(numInt)/denInt
                            err2 = p_hat_mc*(1-p_hat_mc)/denInt
                        elif numInt>denInt:
                            print "numInt>denInt"
                            print "NM1, sample, mlow, mhigh, num, den"
                            print nminus1, mc_Wjets.name, mlow, mhigh, numInt, denInt
                            print
                            p_hat_mc = 1
                            err2 = 0
                        else:
                            if numInt!=0 or denInt!=0:
                                print "ELSE"
                                print "NM1, sample, mlow, mhigh, num, den"
                                print nminus1, mc_Wjets.name, mlow, mhigh, numInt, denInt
                                print
                            p_hat_mc = 0
                            err2 = 0
                        numTot = numTot + numInt*mc_Wjets.partial_weight
#                    print '%s ----- %s ---- %s' % (denTot, denInt, mc.partial_weight)
                        denTot = denTot + denInt*mc_Wjets.partial_weight
                        numsW_Wjets.append(numInt*mc_Wjets.partial_weight)
                        densW_Wjets.append(denInt*mc_Wjets.partial_weight)
                        err2s_Wjets.append(err2)                
                    if denTot == 0:
                	    p_hat = 0
                	    err2sum = 0
                    else:
	                    p_hat = float(numTot)/denTot
                	    err2sum = sum( ((m/denTot)**2 * e2) for m,e2 in zip(densW_Wjets,err2s_Wjets))
                    if err2sum<0:
                        err2sum = 0
                    elif p_hat==0:
                        err2sum = 0
                    p_hats_Wjets.append(p_hat)
                    errsW_Wjets.append(err2sum**0.5)
                    x_Wjets.append(nminus1_num.GetXaxis().GetBinCenter(mbin+1))
                    ex_Wjets.append(nminus1_num.GetXaxis().GetBinWidth(mbin+1)/2)
                eff = ROOT.TGraphAsymmErrors(len(x_Wjets), *[array('d',obj) for obj in (x_Wjets,p_hats_Wjets,ex_Wjets,ex_Wjets,errsW_Wjets,errsW_Wjets)])
            eff.SetTitle(pretty_name)
            ymin, ymax = yrange[nminus1]
            eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin, ymax)
            eff.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
            eff.GetYaxis().SetLabelSize(0.027)
            eff.GetYaxis().SetTitle('n-1 efficiency')
            draw = '2'
            eff.SetLineColor(ROOT.kYellow+2)
            eff.SetFillColor(ROOT.kYellow+2)
            eff.SetFillStyle(3001)
#             lg.AddEntry(eff,'Wjets + Wt + QCD','LF')
            lg.AddEntry(eff,'Single Top','LF')
#             lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
            draw += same
            eff.Draw(draw)
            effs.append(eff)
            same = ' same'
            outfile.cd()
            eff.Write("arp%d"%iarp)
            iarp+=1


#### tt
        print 'tt'
        for entry in to_use_tt[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)

            l = len(mass_range)-1
            nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
            nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))

            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                    if 'data' in entry.name:
                            continue
                #table_wald(entry,nminus1,mass_range)
                    color, fill = styles[entry.name]
                    hnum = entry.histos['NoNo']
                    hden = entry.histos[nminus1]
                    for mbin in range(len(mass_range)):
                            if mbin == (len(mass_range)-1): continue
                            mlow = mass_range[mbin]
                            mhigh = mass_range[mbin+1]
                            num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                            den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                            nminus1_num.SetBinContent(mbin+1, num)
                            nminus1_den.SetBinContent(mbin+1, den)
                    eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
            else:
#             for a,mc in enumerate(entry):
#                 table_wald(mc,nminus1,mass_range)
                p_hats_tt = []
                errsW_tt = []
                x_tt = []
                ex_tt = []
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    numTot = 0
                    denTot = 0
                    err2sum = 0
                    numsW_tt = []
                    densW_tt = []
                    err2s_tt = []
                    for i,mc_tt in enumerate(entry):
                        hnum = mc_tt.histos['NoNo']
                        hden = mc_tt.histos[nminus1]
                        mlow = mass_range[mbin]
                        mhigh = mass_range[mbin+1]
                        numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                        denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                        if numInt<=denInt and denInt>0:
                            p_hat_mc = float(numInt)/denInt
                            err2 = p_hat_mc*(1-p_hat_mc)/denInt
                        elif numInt>denInt:
                            print "numInt>denInt"
                            print "NM1, sample, mlow, mhigh, num, den"
                            print nminus1, mc_tt.name, mlow, mhigh, numInt, denInt
                            print
                            p_hat_mc = 1
                            err2 = 0
                        else:
                            if numInt!=0 or denInt!=0:
                                print "ELSE"
                                print "NM1, sample, mlow, mhigh, num, den"
                                print nminus1, mc_tt.name, mlow, mhigh, numInt, denInt
                                print
                            p_hat_mc = 0
                            err2 = 0
                        numTot = numTot + numInt*mc_tt.partial_weight
#                    print '%s ----- %s ---- %s' % (denTot, denInt, mc.partial_weight)
                        denTot = denTot + denInt*mc_tt.partial_weight
                        numsW_tt.append(numInt*mc_tt.partial_weight)
                        densW_tt.append(denInt*mc_tt.partial_weight)
                        err2s_tt.append(err2)
                    if denTot == 0:
                            p_hat = 0
                            err2sum = 0
                    else:
                            p_hat = float(numTot)/denTot
                            err2sum = sum( ((m/denTot)**2 * e2) for m,e2 in zip(densW_tt,err2s_tt))
                    if err2sum<0:
                        err2sum = 0
                    elif p_hat==0:
                        err2sum = 0
                    p_hats_tt.append(p_hat)
                    errsW_tt.append(err2sum**0.5)
                    x_tt.append(nminus1_num.GetXaxis().GetBinCenter(mbin+1))
                    ex_tt.append(nminus1_num.GetXaxis().GetBinWidth(mbin+1)/2)
                eff = ROOT.TGraphAsymmErrors(len(x_tt), *[array('d',obj) for obj in (x_tt,p_hats_tt,ex_tt,ex_tt,errsW_tt,errsW_tt)])
            eff.SetTitle(pretty_name)
            ymin, ymax = yrange[nminus1]
            eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin, ymax)
            eff.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
            eff.GetYaxis().SetLabelSize(0.027)
            eff.GetYaxis().SetTitle('n-1 efficiency')
            draw = '2'
            eff.SetLineColor(ROOT.kRed+2)
            eff.SetFillColor(ROOT.kRed+2)
            eff.SetFillStyle(3001)
            lg.AddEntry(eff, 'tt','LF')
            #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
            draw += same
            eff.Draw(draw)
            effs.append(eff)
            same = ' same'
            outfile.cd()
            eff.Write("arp%d"%iarp)
            iarp+=1


        lg.Draw()
        name = nminus1+'_split'
        ps.save(name)
        print(nminus1, pretty_name, name)
    # end for name, mass_range in mass_bins:

# ************************************************************************
# Main
# ************************************************************************
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='submits limits for a given mass')
    parser.add_argument('--tag',help='Tagged directory name for input histograms and output plots',default='')
    parser.add_argument('--stats',help='Print stats',default=False)
    parser.add_argument('--lumi',help='Integrated luminosity in pb-1',default='1000')
    parser.add_argument('--do_tight',help='Do tight cuts',default=False)
    args = parser.parse_args()

    # cjsbad - add a cmd-line option for dy only, qcd, zprime, etc
    tag = args.tag
    if not os.path.exists('plots/'):
        raise ValueError('Tagged directory %s does not exist!'%tag)

    draw_mass_test(tag, args.stats, float(args.lumi), args.do_tight)
    
    
print lumi

