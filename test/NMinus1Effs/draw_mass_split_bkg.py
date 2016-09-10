#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

from pprint import pprint
import sys, os
from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)
ROOT.TH1.AddDirectory(0)

from nm1entry import *

print cartella

outfile = ROOT.TFile("whargl_mass.root","recreate")
iarp=0
mcN=0
drawDY = True
do_tight = 'tight' in sys.argv
#psn = 'plots/nminus1effs'
psn = 'plots/TuneP_13_pb_TriggerScale_contributi/' + cartella + '/mass/'
#psn = 'plots/paper2016/mass_onlyDY'
# 'plots' = '/afs/cern.ch/work/c/cschnaib/NMinus1Effs/plots/TAG/'
if do_tight:
    psn += '_tight'
ps = plot_saver(psn, size=(600,600), log=False, pdf=True)
ps.c.SetBottomMargin(0.2)


if do_tight:
    nminus1s = [
        #'TiPt',
        'TiDB',
        'TiGlbChi2',
        'TiIso',
        'TiTkLayers',
        'TiPxHits',
        'TiMuHits',
        'TiMuMatch',
        ]
else:
    nminus1s = [
        'NoPt',
        'NoDB',
        'NoIso',
        'NoTkLayers',
        'NoPxHits',
        'NoMuHits',
        'NoMuMatch',
        'NoVtxProb',
        'NoB2B',
        'NoDptPt',
                #'NoCosm',
        'NoTrgMtch',
        ]


class nm1entry:
    def __init__(self, sample, is_data, lumi):
        if type(sample) == str:
            self.name = sample
            #self.fn = 'data/ana_nminus1_%s_%s.root' % (sample, tipo)
            self.fn = 'data/ana_nminus1_%s.root'%sample
            #self.fn = self.make_fn(sample) if is_data else None
            self.lumi = lumi if is_data else None
            self.is_data = is_data
        else:
            self.name = sample.name
            self.fn = self.make_fn(self.name)
            self.partial_weight = sample.partial_weight
            self.is_data = is_data
        self.prepare_histos()
            
    def make_fn(self, name):
        return 'mc_%s/ana_nminus1_%s.root' % (tipo, name)
    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s + ['NoNo']:
                self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass
                '''
                if nminus1=='NoVtxProb':
                    #self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()
                    self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()
                else:
                    self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass
                '''

    def prepare_histos_sum(self, samples, lumi):
        self.histos = {}
        for nminus1 in nminus1s + ['NoNo']:
            hs = []
            #print '%20s%20s%21s%20s%20s' % ('cut', 'sampe name', 'partial weight', 'scale(ref)','scale(lumi)')
            for sample in samples:
                f = ROOT.TFile(self.make_fn(sample.name))
                if nminus1 == 'NoVtxProb':
                    #h = f.Get(nminus1).Get(cartella).Clone()
                    h = f.Get(nminus1).Get(cartella).Clone()
                else:
                    h = f.Get(nminus1).Get(cartella).Clone()
                #print '%20s%20s%20.15f%20f%20f' % (nminus1, sample.name, sample.partial_weight, refN/refXS, lumiBCD)
                # partial_weight = cross_section * k_factor / Nevents
                if lumi>0:
                    # scale to luminosity for comparision of single dataset to MC
                    h.Scale(sample.partial_weight * lumi) 
                    #print nminus1, sample.name, sample.partial_weight*lumi
                    #print '%20s%20s%20.10f%20f' % (nminus1, sample.name, sample.partial_weight, lumi)
                if lumi<0:
                    # scale to reference cross section/Nevents for comparision of multiple datasets to MC
                    h.Scale(sample.partial_weight * refN / refXS)  
                    #print nminus1, sample.name, sample.partial_weight*refN/refXS
                    #print '%20s%20s%20.10f%20f' % (nminus1, sample.name, sample.partial_weight, refN/refXS)
                hs.append(h)
            hsum = hs[0].Clone()
            for h in hs[1:]:
                hsum.Add(h)
            self.histos[nminus1] = hsum

nolumi = -1
lumi = 13066.191#6293.188#4409.061#4079.344#2231.20#1#2260.881 #2619.44
data = nm1entry('data', True, lumi)#lumiCD )


from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *

raw_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, WZ,ZZ,Wantitop,tW,Wjets,ttbar_lep, dyInclusive50,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200]
use_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, WZ,ZZ,Wantitop,tW,Wjets,ttbar_lep, dyInclusive50,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200]

raw_samples_DY = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000] 
use_samples_DY = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000]

raw_samples_tt = [ttbar_lep]
use_samples_tt = [ttbar_lep]

raw_samples_diboson = [WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, WZ,ZZ]
use_samples_diboson = [WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, WZ,ZZ]

raw_samples_Wjets = [Wantitop,tW, Wjets]
use_samples_Wjets = [Wantitop,tW, Wjets]

mc_samples = [nm1entry(sample,False,lumi) for sample in use_samples]
for mc_sample in mc_samples:
    exec '%s = mc_sample' % mc_sample.name
    
mc_samples_DY = [nm1entry(sample,False,lumi) for sample in use_samples_DY]
for mc_sample_DY in mc_samples_DY:
    exec '%s = mc_sample_DY' % mc_sample_DY.name
    
mc_samples_tt = [nm1entry(sample,False,lumi) for sample in use_samples_tt]
for mc_sample_tt in mc_samples_tt:
    exec '%s = mc_sample_tt' % mc_sample_tt.name

mc_samples_diboson = [nm1entry(sample,False,lumi) for sample in use_samples_diboson]
for mc_sample_diboson in mc_samples_diboson:
    exec '%s = mc_sample_diboson' % mc_sample_diboson.name
    
mc_samples_Wjets = [nm1entry(sample,False,lumi) for sample in use_samples_Wjets]
for mc_sample_Wjets in mc_samples_Wjets:
    exec '%s = mc_sample_Wjets' % mc_sample_Wjets.name

mass_range = [50,120,200,400,800,1400,2300,3500]#,4500, 6000]
to_use = {
    'NoPt':[mc_samples, data],
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
    'NoPt':      (0.5,1.01),
    'NoDB':      (0.5,1.01),
    'NoIso':     (0.1,1.01),
    'NoTkLayers':(0.5,1.01),
    'NoPxHits':  (0.5,1.01),
    'NoMuHits':  (0.5,1.01),
    'NoMuMatch': (0.5,1.01),
    'NoVtxProb': (0.5,1.01),
    'NoB2B':     (0.5,1.01),
    'NoDptPt':   (0.5,1.01),
    'NoTrgMtch': (0.5,1.01),
    }
#global_ymin = 0.
global_ymin = None

def table_wald(entry,nminus1, mass_range):
    print entry.name, nminus1
    hnum = entry.histos['NoNo']
    hden = entry.histos[nminus1]
    #print '%20s%27s%23s%20s%16s%25s%26s' % ('cut', 'mass range','numerator', 'denominator', 'efficiency', '- 68% CL-CP +','68% CL-Wald')
    for mbin in range(len(mass_range)):
        if mbin == (len(mass_range)-1): break
        mlow = mass_range[mbin]
        mhigh = mass_range[mbin+1] 
        num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
        den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
        pcp,lcp,hcp = clopper_pearson(num, den)
        if num==0 and den==0:
            eff = 0
            errw = 0
        else:
            eff = num/den
            if (eff*(1-eff))<0:
                print "what is this"
                print 'NM1, sample, mlow, mhigh, num, den'
                print nminus1, entry.name, mlow, mhigh, num, den
                eff = 1
                errw = 0
            else:
                errw = (eff*(1-eff)/den)**0.5
        if 'data' not in entry.name:
            num = num*entry.partial_weight*lumi
            den = den*entry.partial_weight*lumi
        #print '%20s%15i%15i%20f%20f%15f%15f%15f%23f'     % (nminus1, mlow, mhigh, num, den, eff, eff-lcp, hcp-eff,        errw)
        #print '%20s%15i%15i%20f%20f%15f%15f%15f%15f%16f' % (nminus1, mlow, mhigh, num, den, eff, lcp,     hcp,     eff-errw, eff+errw)
        print ' '
    print '---------------------------------------------'

def table(entry,nminus1, mass_range):
    print entry.name
    hnum = entry.histos['NoNo']
    hden = entry.histos[nminus1]
    print '%20s%27s%23s%20s%20s%22s' % ('cut', 'mass range', 'numerator', 'denominator', 'efficiency', '68% CL')
    for mbin in range(len(mass_range)):
        if mbin == (len(mass_range)-1): break
        mlow = mass_range[mbin]
        mhigh = mass_range[mbin+1] 
        num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
        den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
        e,l,h = clopper_pearson(num, den)
        print '%20s%15i%15i%20f%20f%20f%15f%15f' % (nminus1, mlow, mhigh, num, den, p_hat, p_hat_e, p_hat_e)

ROOT.gStyle.SetTitleX(0.25)
ROOT.gStyle.SetTitleY(0.50)

for nminus1 in nminus1s:
    pretty_name = pretty[nminus1]
    print nminus1, pretty_name
    lg = ROOT.TLegend(0.25, 0.21, 0.91, 0.44)
    lg.SetTextSize(0.03)
    lg.SetFillColor(0)
    lg.SetBorderSize(1)
    
    same = 'A'
    effs = []
    
    
    #TOTALE
    print 'Totale'
    for entry in to_use[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)
       
        l = len(mass_range)-1
        nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
        nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
#         	if 'data' in entry.name:
#         		continue
#         		
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
                    numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                    denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
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
#                    print '%s ----- %s ---- %s' % (denTot, denInt, mc.partial_weight)
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
            #lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (lumi/1000.), 'LP')
            lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (entry.lumi/1000.), 'LP')
        else:
            if len(entry)==9:
                draw = '2'
                eff.SetLineColor(ROOT.kBlue+2)
                eff.SetFillColor(ROOT.kBlue+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'DY MC 80X','LF')
            elif len(entry)==15:
                draw = '2'
                eff.SetLineColor(ROOT.kOrange+2)
                eff.SetFillColor(ROOT.kOrange+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK MC','LF')
            elif len(entry)==1:
                draw = '2'
                eff.SetLineColor(ROOT.kYellow+1)
                eff.SetFillColor(ROOT.kYellow+1)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'Z\'_{#psi} (5 TeV)','LF')
            elif len(entry)==10:
                draw = '2'
                eff.SetLineColor(ROOT.kViolet)
                eff.SetFillColor(ROOT.kViolet)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'QCD MC','LF')
            elif len(entry)==11:
                draw = '2'
                eff.SetLineColor(ROOT.kBlue+2)
                eff.SetFillColor(ROOT.kBlue+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'QCD + W+jets MC','LF')
            elif len(entry)==7:
                draw = '2'
                eff.SetLineColor(ROOT.kRed+2)
                eff.SetFillColor(ROOT.kRed+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK no DY MC','LF')
            elif len(entry)==6:
                draw = '2'
                eff.SetLineColor(ROOT.kRed+2)
                eff.SetFillColor(ROOT.kRed+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK','LF')
            elif len(entry)==len(use_samples):
                draw = '2'
                eff.SetLineColor(ROOT.kGreen+2)
                eff.SetFillColor(ROOT.kGreen+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'All Simulation','LF')
            else:
                draw = '2'
                eff.SetLineColor(ROOT.kBlue+2)
                eff.SetFillColor(ROOT.kBlue+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK+QCD MC','LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        outfile.cd()
        eff.Write("arp%d"%iarp)
        iarp+=1
    # end for entry in to_use[name]: # entry is a specific sample
    
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
        eff.SetLineColor(7)
        eff.SetFillColor(7)
        eff.SetFillStyle(1001)
        lg.AddEntry(eff,'DY','LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        outfile.cd()
        eff.Write("arp%d"%iarp)
        iarp+=1

    #DIBOSON
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
        eff.SetFillStyle(1001)
        lg.AddEntry(eff,'Diboson','LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        outfile.cd()
        eff.Write("arp%d"%iarp)
        iarp+=1

    #tt
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
        eff.SetFillStyle(1001)
        lg.AddEntry(eff, tt,'LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        outfile.cd()
        eff.Write("arp%d"%iarp)
        iarp+=1

    #Wjets
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
                print p_hat
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
        eff.SetFillStyle(1001)
        lg.AddEntry(eff,'Wjets + Wt','LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        outfile.cd()
        eff.Write("arp%d"%iarp)
        iarp+=1

        
    t = ROOT.TPaveLabel(0.50, 0.525, 0.90, 0.625, categoria, 'brNDC')
    t.SetTextFont(42)
    t.SetTextSize(0.5)
    t.SetBorderSize(0)
    t.SetFillColor(0)
    t.SetFillStyle(0)
    tt = ROOT.TPaveLabel(0.50, 0.425, 0.90, 0.525, categoria_1, 'brNDC')
    tt.SetTextFont(42)
    tt.SetTextSize(0.5)
    tt.SetBorderSize(0)
    tt.SetFillColor(0)
    tt.SetFillStyle(0)
    tt.Draw() 
    t.Draw() 
    
    lg.Draw()
    ps.save(nminus1+'_mass')
    print
print cartella, tipo
# end for name, mass_range in mass_bins:
