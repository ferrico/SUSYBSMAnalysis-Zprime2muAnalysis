#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

from pprint import pprint
import sys, os
from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *

from nm1entry import *

set_zp2mu_style()
#ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
#ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)
ROOT.TH1.AddDirectory(0)


outfile = ROOT.TFile("whargl_mass_stack.root","recreate")
iarp=0
do_tight = 'tight' in sys.argv
#psn = 'plots/nminus1effs'
psn = 'plots/TuneP_13_pb_TriggerScale_contributi/' + cartella + '/mass_stack/'
# 'plots' = '/afs/cern.ch/work/c/cschnaib/NMinus1Effs/plots/TAG/'
if do_tight:
    psn += '_tight'
ps = plot_saver(psn, size=(800,600),log=True, pdf_log=True, name='')
ps.c.SetTopMargin(0.075)
ps.c.SetBottomMargin(0.1)
ps.c.SetRightMargin(0.05)
#ps.c.SetLeftMargin(0.03)
ps.c.SetLogy(1)


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

pretty = {
    'NoPt': 'p_{T} > 53 GeV',
    'NoTkLayers': '# tk lay > 5',
    'NoPxHits': '# px hits > 0',
    'NoMuStns': '# mu segs > 1',
    'NoDB': '|dxy| < 0.2',
    'NoGlbChi2': 'glb #chi^{2}/ndf < 10',
    'NoTkMuon': 'isTrackerMuon',
    'NoMuHits': '# mu hits > 0',
    'NoMuMatch': '# matched stations > 1',
    'NoCosm': 'anti-cosmic',
    'NoTrgMtch': 'HLT match',
    'NoB2B': 'back-to-back',
    'NoVtxProb': '#chi^{2} #mu#mu vtx < 20',
    'NoDptPt': 'dpT/pT',
    'NoIso': 'rel. tk. iso.',
    #'data': 'Data, %.1f fb^{-1}',
    'data': 'Data, %.1f fb^{-1}, MuonOnly',
    'dataB': 'Data RunB, %.1f fb^{-1}, MuonOnly',
    'dataCD': 'Data RunC+D, %.1f fb^{-1}, MuonOnly',
    'dataBCD': 'Data RunB+C+D, %.1f fb^{-1}, MuonOnly',
    'mcsum_lumi': 'Simulation',
    'mcsum_ref': 'Simulation',
    'mc50m120_lumi': 'Simulation 60 < M < 120 GeV',
    'mc50m120_ref': 'Simulation 60 < M < 120 GeV',
    'mc120m_lumi': 'Simulation M > 120 GeV',
    'mc120m_ref': 'Simulation M > 120 GeV',
    'mc800m2300_lumi': 'Simulation 800 < M < 2300 GeV',
    'mc800m2300_ref': 'Simulation 800 < M < 2300 GeV',
    'mc400m2300_lumi': 'Simulation 400 < M < 2300 GeV',
    'mc400m2300_ref': 'Simulation 400 < M < 2300 GeV',
    'zmumu': 'Z#rightarrow#mu#mu, 60 < M < 120 GeV',
    'dy120_c1': 'DY#rightarrow#mu#mu, M > 120 GeV',
    'dy200_c1': 'DY#rightarrow#mu#mu, M > 200 GeV',
    'dy500_c1': 'DY#rightarrow#mu#mu, M > 500 GeV',
    'dy1000_c1': 'DY#rightarrow#mu#mu, M > 1000 GeV',
    'dy50': 'DY#rightarrow#mu#mu madgraph',
#    'dy50': 'DY#rightarrow#mu#mu, M > 50 GeV',
    'dy50to120': 'DY#rightarrow#mu#mu 50 < m < 120 GeV',
    'dy120to200': 'DY#rightarrow#mu#mu 120 < m < 200 GeV',
    'dy200to400': 'DY#rightarrow#mu#mu 200 < m < 400 GeV',
    'dy400to800': 'DY#rightarrow#mu#mu 400 < m < 800 GeV',
    'dy800to1400': 'DY#rightarrow#mu#mu 800 < m < 1400 GeV',
    'dy1400to2300': 'DY#rightarrow#mu#mu 1400 < m < 2300 GeV',
    'dy2300to3500': 'DY#rightarrow#mu#mu 2300 < m < 3500 GeV',
    'dy3500to4500': 'DY#rightarrow#mu#mu 3500 < m < 4500 GeV',
    'dy4500to6000': 'DY#rightarrow#mu#mu 4500 < m < 6000 GeV',
    'dy50_startup': 'DY#rightarrow#mu#mu startup',
    'ttbar': 't#bar{t}',
    'ttbar_pow': 't#bar{t} powheg',
    'ttbar_startup': 't#bar{t} startup',
    'WWinclusive': 'WW: 50 < m < 200',
    'WW200to600': 'WW: 200 < m < 600',
    'WW600to1200': 'WW: 600 < m < 1200',
    'WW1200to2500': 'WW: 1200 < m < 2500',
    'WW2500': 'WW: m > 2500',
    'ZZ': 'ZZ',
    'WZ' : 'WZ',
    'dyInclusive' : 'TauTau',
    'inclmu15': 'QCD',
    'zssm1000': 'Z\' SSM, M=1000 GeV',
    'zpsi5000': 'Z\'_{#psi}, M=5000 GeV',
    'zpsi5000_m1TeV': 'Z\'_{#psi}, M=5000 GeV',
    'zpsi5000_1m3TeV': 'Z\'_{#psi}, M=5000 GeV',
    'zpsi5000_3mTeV': 'Z\'_{#psi}, M=5000 GeV',
    '60m120_BCD': '60 < m < 120 GeV',
    '60m120_CD': '60 < m < 120 GeV',
    '60m120': '60 < m < 120 GeV',
    '70m110': '70 < m < 110 GeV',
    '120m200': '120 < m < 200 GeV', 
    '200m400': '200 < m < 400 GeV',
    '400m600': '400 < m < 600 GeV',
    '200m': 'm > 200 GeV',
    '50m': 'm > 50 GeV',
    '70m': 'm > 70 GeV',
    '120m_BCD': 'm > 120 GeV',
    '120m_CD': 'm > 120 GeV',
    '120m': 'm > 120 GeV',
    'DY120to200powheg': 'DY#rightarrow#mu#mu 120 < m < 200 GeV',
    'DY200to400powheg': 'DY#rightarrow#mu#mu 200 < m < 400 GeV',
    'DY400to800powheg': 'DY#rightarrow#mu#mu 400 < m < 800 GeV',
    'DY800to1400powheg': 'DY#rightarrow#mu#mu 800 < m < 1400 GeV',
    'dy1400to2300': 'DY#rightarrow#mu#mu 1400 < m < 2300 GeV',
    '400m800' : '400 < m < 800 GeV',
    '800m1400': '800 < m < 1400 GeV',
    '1400m2300':'1400 < m < 2300 GeV',
    '800m2300':'800 < m < 2300 GeV',
    '400m2300':'400 < m < 2300 GeV',
    'all_lumi':'Simulation M > 120 GeV',
    'all_ref':'Simulation M > 120 GeV',
    '120m1400':'120 < M < 1400 GeV',
    }

class nm1entry:
    def __init__(self, sample, is_data, lumi):
        if type(sample) == str:
            self.name = sample
            #self.fn = 'data/ana_nminus1_%s_%s.root' % (sample, tipo)
            self.fn = 'data/ana_nminus1_%s.root'%sample
            #self.fn = self.make_fn(sample,is_data) if is_data else None
            self.lumi = lumi if is_data else None
            self.is_data = is_data
        else:
            self.name = sample.name
            self.fn = self.make_fn(self.name,is_data)
            self.partial_weight = sample.partial_weight
            self.is_data = is_data
        self.prepare_histos()
            
    def make_fn(self, name, is_data):
        '''
        if is_data==True:
            return 'data/ana_nminus1_%s.root' %name
        else:
        '''
        return 'mc_%s/ana_nminus1_%s.root' % (tipo, name)
    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s + ['NoNo']:
                if nminus1=='NoVtxProb':
                    #self.histos[nminus1] = f.Get(nminus1).Get('DileptonMass').Clone()
                    self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass
                else:
                    self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass

    def prepare_histos_sum(self, samples, lumi):
        self.histos = {}
        for nminus1 in nminus1s + ['NoNo']:
            #hs = []
            #print '%20s%20s%21s%20s%20s' % ('cut', 'sampe name', 'partial weight', 'scale(ref)','scale(lumi)')
            for sample in samples:
                f = ROOT.TFile(self.make_fn(sample.name,is_data))
                if nminus1 == 'NoVtxProb':
                    #h = f.Get(nminus1).Get('DileptonMass').Clone()
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
                #hs.append(h)
            #hsum = hs[0].Clone()
            #for h in hs[1:]:
            #    hsum.Add(h)
            self.histos[nminus1] = h


nolumi = -1
lumi = 13066.191#6293.188#4409.061#4079.344#2231.20#1#2260.881 #2619.44
data = nm1entry('data', True, lumi)#lumiCD )

from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
raw_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, WZ,ZZ,Wantitop,tW,Wjets,ttbar_lep, dyInclusive50,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200]#,
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

#yrange = {
#   'sample':    (ymin,ymax),
#    'NoPt':      (0.00,1.01),
#    'NoDB':      (0.95,1.001),
#    'NoIso':     (0.6,1.01),
#    'NoTkLayers':(0.95,1.001),
#    'NoPxHits':  (0.95,1.001),
#    'NoMuHits':  (0.95,1.001),
#    'NoMuMatch': (0.35,1.005),
#    'NoVtxProb': (0.90,1.001),
#    'NoB2B':     (0.95,1.001),
#    'NoDptPt':   (0.95,1.001),
#    'NoTrgMtch': (0.95,1.001),
#    }
#global_ymin = 0.
global_ymin = None

def table_wald(entry,nminus1, mass_range):
    print entry.name
    hnum = entry.histos['NoNo']
    hden = entry.histos[nminus1]
    print '%20s%27s%23s%20s%16s%25s%26s' % ('cut', 'mass range','numerator', 'denominator', 'efficiency', '- 68% CL-CP +','68% CL-Wald')
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
                print nminus1, entry.name, mlow, mhigh, num, den
            else:
                errw = (eff*(1-eff)/den)**0.5
        print '%20s%15i%15i%20f%20f%15f%15f%15f%23f'     % (nminus1, mlow, mhigh, num, den, eff, eff-lcp, hcp-eff,        errw)
        print '%20s%15i%15i%20f%20f%15f%15f%15f%15f%16f' % (nminus1, mlow, mhigh, num, den, eff, lcp,     hcp,     eff-errw, eff+errw)
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

#ROOT.gStyle.SetTitleX(0.25)
#ROOT.gStyle.SetTitleY(0.50)

for nminus1 in nminus1s:
    pretty_name = pretty[nminus1]
    print nminus1, pretty_name
    lg = ROOT.TLegend(0.45, 0.55, 0.91, 0.9)
    lg.SetTextSize(0.03)
    lg.SetFillColor(0)
    lg.SetBorderSize(1)
    
    same = 'A'
    effs = []

    stack_num = ROOT.THStack('hs','')
    stack_den = ROOT.THStack('hs','')



    for entry in to_use[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)

        l = len(mass_range)-1

        data_num = ROOT.TH1F('num', '', l, array('f',mass_range))
        data_den = ROOT.TH1F('den', '', l, array('f',mass_range))

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
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
                #nminus1_num.SetBinContent(mbin+1, num)
                data_num.SetBinContent(mbin+1, num)
                data_den.SetBinContent(mbin+1, den)
            #eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)

            draw = 'P'
            data_num.SetLineColor(color)
            data_num.SetMarkerStyle(20)
            data_num.SetMarkerSize(1.05)
            data_num.SetMarkerColor(color)
            data_den.SetLineColor(color)
            data_den.SetMarkerStyle(20)
            data_den.SetMarkerSize(1.05)
            data_den.SetMarkerColor(color)
            lg.AddEntry(data_den, pretty.get(entry.name, entry.name) % (entry.lumi/1000.), 'LP')



    for entry in to_use_DY[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)

        l = len(mass_range)-1

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
        	continue;
        else:
            #for a,mc in enumerate(entry):
             #   table_wald(mc,nminus1,mass_range)
            for i,mc in enumerate(entry):
                nminus1_den_MC = ROOT.TH1F('den', '', l, array('f',mass_range))
                nminus1_num_MC = ROOT.TH1F('num', '', l, array('f',mass_range))
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    hnum = mc.histos['NoNo']
                    hden = mc.histos[nminus1]
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                    #den*mc.partial_weight*lumiBCD
                    nminus1_num_MC.SetBinContent(mbin+1, num)
                    nminus1_den_MC.SetBinContent(mbin+1, den)
                color, fill = styles[mc.name]
                nminus1_num_MC.SetLineColor(7)
                nminus1_num_MC.SetFillColor(7)
                nminus1_num_MC.SetMinimum(1E-7)
                #nminus1_num_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_num_MC.Scale(mc.partial_weight*lumi)
                stack_num.Add(nminus1_num_MC)
                nminus1_den_MC.SetLineColor(7)
                nminus1_den_MC.SetFillColor(7)
                #nminus1_den_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_den_MC.Scale(mc.partial_weight*lumi)
                stack_den.Add(nminus1_den_MC)
                #lg.AddEntry(nminus1_den_MC,pretty.get(mc.name,mc.name),"F")
            lg.AddEntry(nminus1_den_MC,"Drell-Yan",'F' )
            
            
    for entry in to_use_diboson[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)

        l = len(mass_range)-1

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
        	continue;
        else:
            #for a,mc in enumerate(entry):
             #   table_wald(mc,nminus1,mass_range)
            for i,mc in enumerate(entry):
                nminus1_den_MC = ROOT.TH1F('den', '', l, array('f',mass_range))
                nminus1_num_MC = ROOT.TH1F('num', '', l, array('f',mass_range))
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    hnum = mc.histos['NoNo']
                    hden = mc.histos[nminus1]
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                    #den*mc.partial_weight*lumiBCD
                    nminus1_num_MC.SetBinContent(mbin+1, num)
                    nminus1_den_MC.SetBinContent(mbin+1, den)
                color, fill = styles[mc.name]
                nminus1_num_MC.SetLineColor(ROOT.kBlue+2)
                nminus1_num_MC.SetFillColor(ROOT.kBlue+2)
                nminus1_num_MC.SetMinimum(1E-7)
                #nminus1_num_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_num_MC.Scale(mc.partial_weight*lumi)
                stack_num.Add(nminus1_num_MC)
                nminus1_den_MC.SetLineColor(ROOT.kBlue+2)
                nminus1_den_MC.SetFillColor(ROOT.kBlue+2)
                #nminus1_den_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_den_MC.Scale(mc.partial_weight*lumi)
                stack_den.Add(nminus1_den_MC)
                #lg.AddEntry(nminus1_den_MC,pretty.get(mc.name,mc.name),"F")
            lg.AddEntry(nminus1_den_MC,"Diboson",'F' )
                


    for entry in to_use_tt[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)

        l = len(mass_range)-1

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
        	continue;
        else:
            #for a,mc in enumerate(entry):
             #   table_wald(mc,nminus1,mass_range)
            for i,mc in enumerate(entry):
                nminus1_den_MC = ROOT.TH1F('den', '', l, array('f',mass_range))
                nminus1_num_MC = ROOT.TH1F('num', '', l, array('f',mass_range))
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    hnum = mc.histos['NoNo']
                    hden = mc.histos[nminus1]
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                    #den*mc.partial_weight*lumiBCD
                    nminus1_num_MC.SetBinContent(mbin+1, num)
                    nminus1_den_MC.SetBinContent(mbin+1, den)
                color, fill = styles[mc.name]
                nminus1_num_MC.SetLineColor(ROOT.kRed+2)
                nminus1_num_MC.SetFillColor(ROOT.kRed+2)
                nminus1_num_MC.SetMinimum(1E-7)
                #nminus1_num_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_num_MC.Scale(mc.partial_weight*lumi)
                stack_num.Add(nminus1_num_MC)
                nminus1_den_MC.SetLineColor(ROOT.kRed+2)
                nminus1_den_MC.SetFillColor(ROOT.kRed+2)
                #nminus1_den_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_den_MC.Scale(mc.partial_weight*lumi)
                stack_den.Add(nminus1_den_MC)
                #lg.AddEntry(nminus1_den_MC,pretty.get(mc.name,mc.name),"F")
            lg.AddEntry(nminus1_den_MC, tt,'F' )
  
            
    for entry in to_use_Wjets[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)

        l = len(mass_range)-1

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
        	continue;
        else:
            #for a,mc in enumerate(entry):
             #   table_wald(mc,nminus1,mass_range)
            for i,mc in enumerate(entry):
                nminus1_den_MC = ROOT.TH1F('den', '', l, array('f',mass_range))
                nminus1_num_MC = ROOT.TH1F('num', '', l, array('f',mass_range))
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    hnum = mc.histos['NoNo']
                    hden = mc.histos[nminus1]
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                    #den*mc.partial_weight*lumiBCD
                    nminus1_num_MC.SetBinContent(mbin+1, num)
                    nminus1_den_MC.SetBinContent(mbin+1, den)
                color, fill = styles[mc.name]
                nminus1_num_MC.SetLineColor(ROOT.kYellow+2)
                nminus1_num_MC.SetFillColor(ROOT.kYellow+2)
                nminus1_num_MC.SetMinimum(1E-7)
                #nminus1_num_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_num_MC.Scale(mc.partial_weight*lumi)
                stack_num.Add(nminus1_num_MC)
                nminus1_den_MC.SetLineColor(ROOT.kYellow+2)
                nminus1_den_MC.SetFillColor(ROOT.kYellow+2)
                #nminus1_den_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_den_MC.Scale(mc.partial_weight*lumi)
                stack_den.Add(nminus1_den_MC)
                #lg.AddEntry(nminus1_den_MC,pretty.get(mc.name,mc.name),"F")
            lg.AddEntry(nminus1_den_MC,"Wjets + Wt",'F' )





    #stack_num.SetMinimum(1E-7)
    t = ROOT.TPaveLabel(0.50, 0.525, 0.90, 0.625, categoria, 'brNDC')
    t.SetTextFont(42)
    t.SetTextSize(0.5)
    t.SetBorderSize(0)
    t.SetFillColor(0)
    t.SetFillStyle(0)
    tt = ROOT.TPaveLabel(0.50, 0.425, 0.90, 0.525, categoria, 'brNDC')
    tt.SetTextFont(42)
    tt.SetTextSize(0.5)
    tt.SetBorderSize(0)
    tt.SetFillColor(0)
    tt.SetFillStyle(0)

    stack_num.Draw("hist")
    data_num.Draw("pe1same")
    tt.Draw()
    t.Draw() 
    lg.Draw()
    outfile.cd()
    #stack_num.SetMinimum(0.1)
    stack_num.SetTitle("All Selection Applied")
    stack_num.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    #stack_num.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
    stack_num.GetYaxis().SetTitle("Events")
    ps.save(nminus1+'_stack_mass_num')
    print
    
    stack_den.Draw("hist")
    data_den.Draw("pe1same")
    tt.Draw()
    t.Draw() 
    lg.Draw()
    outfile.cd()
    stack_den.SetMinimum(0.1)
    stack_den.SetTitle("N - ("+pretty_name+")")
    stack_den.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    #stack_den.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
    stack_den.GetYaxis().SetTitle("Events")
    ps.save(nminus1+'_stack_mass_den')
    print
    
print cartella, tipo
# end for name, mass_range in mass_bins:
