#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

from pprint import pprint
import sys, os
from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
#ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
#ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)
ROOT.TH1.AddDirectory(0)

from nm1entry import *

print cartella

outfile = ROOT.TFile("whargl_mass_stack_2.root","recreate")
iarp=0
do_tight = 'tight' in sys.argv
#psn = 'plots/nminus1effs'
# psn = 'plots/' + period + categoria_1 + cartella + '/mass_stack/'
psn = 'plots/perGiovanna_Nminus1/mass_stack/' + cartella 
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


class nm1entry:
    def __init__(self, sample, is_data, lumi):
        if type(sample) == str:
            self.name = sample
            self.fn = '/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/NMinus1Effs/data/Run2016MuonsOnly/ana_nminus1_%s.root' %sample
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
            return 'mc/ana_nminus1_%s.root' %name
        else:
        '''
        return 'mc/ana_nminus1_%s.root' % name

    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s + ['NoNo']:
                self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass
                '''
                if nminus1=='NoVtxProb':
                    #self.histos[nminus1] = f.Get(nminus1).Get('DileptonMass').Clone()
                    self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()
                else:
                    self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass
                '''

#     def prepare_histos(self):
#         self.histos = {}
#         if self.fn is not None:
#             f = ROOT.TFile(self.fn)
#             for nminus1 in nminus1s + ['NoNo']:
# #                 if nminus1=='NoVtxProb':
#                     #self.histos[nminus1] = f.Get(nminus1).Get('DileptonMass').Clone()
#                     self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass
# #                 else:
# #                     self.histos[nminus1] = f.Get(nminus1).Get(cartella).Clone()#DileptonMass

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

#data, lumi = nm1entry('data', True), 242.8 # lumi in pb
nolumi = -1
#lumiB = 50.7 
#lumi = 13981.143#13066.191#6293.188#4079.344#2231.20#1#2260.881 #2619.44

#lumi = 36781.461
#lumi = 13081.143

#lumiD = 2572.19
#lumiBCD = 2660.14
#lumiBCD = 2800.
#dataB = nm1entry('dataB', True, lumiB)#lumiB )
data = nm1entry('data', True, lumi)#lumiCD )
#dataBCD = nm1entry('dataBCD', True, lumiBCD)#lumiCD )
#dataD = nm1entry('dataD', True, lumiD )
#mcsum_lumi = nm1entry('mcsum_lumi',False,lumiBCD)
#mcsum_ref = nm1entry('mcsum_ref',False,nolumi)
#DYmc = nm1entry('DYmc',False,lumiBCD)
#nonDYmc = nm1entry('nonDYmc',False,lumiBCD)
#Wjets = nm1entry('Wjets',False,lumiBCD)

from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *			
use_samples = [ Wantitop, tW, 
				ttbar_lep50to500, ttbar_lep_500to800, ttbar_lep_800to1200, ttbar_lep_1200to1800, ttbar_lep_1800,
				dyInclusive50,
				Wjets, 
				qcd80to120, qcd120to170, qcd170to300, qcd300to470, qcd470to600, qcd600to800, 
				qcd800to1000, qcd1000to1400, qcd1400to1800, qcd1800to2400, qcd2400to3200, qcd3200,
				WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, 
				ZZ, WZ, 
				ZZ_ext, WZ_ext, 
				dy50to120, 
				dy120to200, dy200to400, dy400to800, dy800to1400, 
				dy1400to2300, dy2300to3500, dy3500to4500, dy4500to6000
				]


refXS = dy50to120.cross_section
refN = dy50to120.nevents
print lumi, refN/refXS

# All MC samples
# lumi
#mc_samples.prepare_histos_sum(use_samples,lumiBCD)
mc_samples = [nm1entry(sample,False,lumi) for sample in use_samples]
for mc_sample in mc_samples:
    exec '%s = mc_sample' % mc_sample.name
# ref
#mcsum_ref.prepare_histos_sum(use_samples, nolumi) 
#mc_samples_ref = [nm1entry(sample,False,nolumi) for sample in use_samples]
#for mc_sample in mc_samples_ref:
#    exec '%s = mc_sample' % mc_sample.name

#bin_width = 20
#maxX = 2500
#minX = 60
#nBins = (maxX-minX)/bin_width
#mass_range = []
#for i in range(3,nBins):
#    ibin = i*bin_width
#    mass_range.append(ibin)
#print mass_range

# mass_range = [40,60,80,100,120,150,180,240,320,500,1000,2000,2500]
# mass_range = [50,120,200,300,400,500,600,800,1400,2300]#,3500]#,4500, 6000]

# mass_range = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
# 			11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
# 			21, 22, 23, 24, 25, 26, 27, 28, 29, 30] #DimuonMassVtx_chi2
# mass_range = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 
# 				220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 
# 				420, 440, 460, 480, 500, 520, 540, 560, 580, 600,
# 				620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 
# 				820, 840, 860, 880, 900, 920, 940, 960, 980, 1000] # DileptonPt
# mass_range = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.087, 0.09, 0.095, 0.1, 
# 			  0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 
# 			  0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.3] #RelIsoSumPt
# mass_range = [0,  0.02, 0.04,  0.06,  0.08, 0.1, 
# 			   0.12, 0.14, 0.16, 0.18, 0.2, 
# 			   0.22, 0.24, 0.26, 0.28, 0.3,
# 			   0.32, 0.34, 0.36, 0.38, 0.4, 
# 			   0.42, 0.44, 0.46, 0.48, 0.5,
# 			   0.52, 0.54, 0.56, 0.58]#, 0.6, 
# 			   0.62, 0.64, 0.66, 0.68, 0.7,
#			   0.72, 0.74, 0.76, 0.78, 0.8,
#			   0.82, 0.84, 0.86, 0.88, 0.9, 
#			   0.92, 0.94, 0.96, 0.98, 1, 1.02] #RelIsoSumPt
mass_range = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
				10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
				20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
				30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
			40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
			50]#, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60] #NoVtx
# mass_range = [40, 60, 80, 100, 120, 140, 160, 180, 200,400,800,1400]
# mass_range = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
# 				0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1] #DimuonMuonPtErrOverPt
# mass_range = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 			  11, 12, 13, 14, 15, 16, 17, 18, 19, 20] #NTkLayers
# mass_range = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #NPxHits
# mass_range = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]#,
# 			  11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
# 			  21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
# 			  31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
# 			  41, 42, 43, 44, 45, 46, 47, 48, 49, 50] #NVertices
#50,120,200,400,800,1400,2300,3500]#,4500, 6000] 

to_use = {

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
#     lg = ROOT.TLegend(0.35, 0.15, 0.71, 0.5)  #vtx chi2, pixelHits, 
#     lg = ROOT.TLegend(0.55, 0.55, 0.91, 0.9)  #sgimaPtSuPt, 
    lg = ROOT.TLegend(0.45, 0.15, 0.81, 0.45)  #tracker
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
        
        NUM = 0
        DEN = 0
        

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
#                 nminus1_num.SetBinContent(mbin+1, num)
                data_num.SetBinContent(mbin+1, num)
                data_den.SetBinContent(mbin+1, den)
            #eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
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
                    NUM = NUM + num
                    DEN = DEN + den
#                     print mc.name, mbin, num, NUM, den, DEN
                    nminus1_num_MC.SetBinContent(mbin+1, num)
                    nminus1_den_MC.SetBinContent(mbin+1, den)
#                 print mc.name, NUM, DEN
                color, fill = styles[mc.name]
                nminus1_num_MC.SetLineColor(color)
                nminus1_num_MC.SetFillColor(color)
                nminus1_num_MC.SetMinimum(1E-7)
                #nminus1_num_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_num_MC.Scale(mc.partial_weight*lumi)
                stack_num.Add(nminus1_num_MC)
                nminus1_den_MC.SetLineColor(color)
                nminus1_den_MC.SetFillColor(color)
                #nminus1_den_MC.Scale(mc.partial_weight*refN/refXS)
                nminus1_den_MC.Scale(mc.partial_weight*lumi)
                stack_den.Add(nminus1_den_MC)
#                 lg.AddEntry(nminus1_den_MC,pretty.get(mc.name,mc.name),"F")
                if mc.name == 'dy50to120':
                    lg.AddEntry(nminus1_den_MC,"#gamma/Z #rightarrow #mu^{+}#mu^{-}",'F' )
                elif mc.name == 'dyInclusive50':
                    lg.AddEntry(nminus1_den_MC,"#gamma/Z #rightarrow #tau^{+}#tau^{-}",'F' )
                elif 'ttbar_lep50to500' in mc.name:
                    lg.AddEntry(nminus1_den_MC,"t#bar{t}","F")
#                elif mc.name == 'WWinclusive':
                elif mc.name == 'WZ':
                    lg.AddEntry(nminus1_den_MC,"WZ","F")
                elif mc.name == 'ZZ':
                    lg.AddEntry(nminus1_den_MC,"ZZ","F")
                elif 'WWinclusive' in mc.name:
                    lg.AddEntry(nminus1_den_MC,"WW","F")
                elif mc.name == 'Wjets':
                    lg.AddEntry(nminus1_den_MC,"W+jets","F")
                elif 'tW' == mc.name:
                    lg.AddEntry(nminus1_den_MC,"Single Top","F")
                elif 'qcd80to120' in mc.name:
                    lg.AddEntry(nminus1_den_MC,"QCD","F")

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
            outfile.cd()
       	    data_num.Write(nminus1+'_data_num')#(nminus1+'_stack_mass_num')
       	    data_den.Write(nminus1+'_data_den')#(nminus1+'_stack_mass_den')
            draw = 'P'
            data_num.SetLineColor(color)
            data_num.SetMarkerStyle(20)
            data_num.SetMarkerSize(1.05)
            data_num.SetMarkerColor(color)
            data_den.SetLineColor(color)
            data_den.SetMarkerStyle(20)
            data_den.SetMarkerSize(1.05)
            data_den.SetMarkerColor(color)
#             lg.AddEntry(data, pretty.get(entry.name, entry.name) % (lumi/1000.), 'LP')
            lg.AddEntry(data_den, pretty.get(entry.name, entry.name) % (entry.lumi/1000.), 'LP')
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
#     t.Draw()
#     tt.Draw()
    lg.Draw()
    outfile.cd()
    stack_num.SetMinimum(0.00005)
    stack_num.SetTitle("All Selection Applied")
    stack_num.GetXaxis().SetTitle(axisX)
    #stack_num.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
    stack_num.GetYaxis().SetTitle("Events")
    stack_num.Write(nminus1+'_stack_mass_num')
    ps.save(nminus1+'_stack_mass_num')
    print
    
    stack_den.Draw("hist")
    data_den.Draw("pe1same")
#     t.Draw()
#     tt.Draw()
    lg.Draw()
    outfile.cd()
    stack_den.SetMinimum(0.00005)
    stack_den.SetTitle("N - ("+pretty_name+")")
    stack_den.GetXaxis().SetTitle(axisX)
    #stack_den.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
    stack_den.GetYaxis().SetTitle("Events")
    stack_den.Write(nminus1+'_stack_mass_den')
    ps.save(nminus1+'_stack_mass_den')
    print
    
print cartella, axisX, lumi
# end for name, mass_range in mass_bins:
