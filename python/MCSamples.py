#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn, files_from_dbs
from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import dataset_from_publish_log

class sample(object):
    def __init__(self, name, nice_name, dataset, nevents, color, syst_frac, cross_section, k_factor=1, filenames=None, scheduler='condor', hlt_process_name='HLT', ana_dataset=None, is_madgraph=False, is_zprime=False):
        self.name = name
        self.nice_name = nice_name
        self.dataset = dataset
        self.nevents = nevents
        self.color = color
        self.syst_frac = syst_frac
        self.cross_section = cross_section
        self.k_factor = k_factor
        self.filenames_ = filenames
        self.scheduler = scheduler
        self.hlt_process_name = hlt_process_name
        self.ana_dataset = ana_dataset
        self.is_madgraph = is_madgraph
        self.is_zprime = is_zprime

    @property
    def partial_weight(self):
        return self.cross_section / float(self.nevents) * self.k_factor # the total weight is partial_weight * integrated_luminosity

    @property
    def filenames(self):
        # Return a list of filenames for running the histogrammer not
        # using crab.
        if self.filenames_ is not None:
            return self.filenames_
        return files_from_dbs(self.ana_dataset, ana02=True)

    def __getitem__(self, key):
        return getattr(self, key)

    def _dump(self, redump_existing=False):
        dst = os.path.join('/uscmst1b_scratch/lpc1/3DayLifetime/tucker', self.name) 
        os.system('mkdir ' + dst)
        for fn in self.filenames:
            print fn
            if redump_existing or not os.path.isfile(os.path.join(dst, os.path.basename(fn))):
                os.system('dccp ~%s %s/' % (fn,dst))

class tupleonlysample(sample):
    def __init__(self, name, dataset, scheduler='condor', hlt_process_name='HLT'):
        super(tupleonlysample, self).__init__(name, 'dummy', dataset, 1, 1, 1, 1, scheduler=scheduler, hlt_process_name=hlt_process_name)

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV for xsecs (all below in pb)
# Single-top cross sections are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma
# K factor for Drell-Yan samples is the ratio of the NNLO to POWHEG cross sections for M > 20 GeV bin, 1915/1871=1.024
samples = [
		#sample('dyInclusive50',          'DYInclusive50', '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 28696958, 209 , 1., 6025.2,    k_factor=1., is_madgraph=True),  #(I will update the number of events)
		
		sample('dy50to120',          'DY50to120', '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 2967200, 209 , 1., 1975,   k_factor=1.006),#NLO xs and k-factor applied to reach NLO
        sample('dy120to200',     'DY120to200', '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 99200, 210, 1., 19.32, k_factor=1.006),#mcm 19.32
        sample('dy200to400',  'DY200to400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 211, 1., 2.731, k_factor=1.006),#mcm 2.731
        sample('dy400to800',  'DY400to800', '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/AODSIM', 100000, 212, 1., 0.241, k_factor=1.006),
        sample('dy800to1400',     'DY800to1400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 72, 1., 0.01678, k_factor=1.006),
        sample('dy1400to2300',          'DY1400to2300', '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 70 , 1., 0.00139,    k_factor=1.006),
        sample('dy2300to3500',          'DY2300to3500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 99200, 70 , 1., 0.00008948,    k_factor=1.006),
        sample('dy3500to4500',          'DY3500to4500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 70 , 1., 0.0000041,    k_factor=1.006),
        sample('dy4500to6000',          'DY4500to6000', '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 70 , 1., 4.56E-7,    k_factor=1.006),

        sample('WZ',        'WZ', '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 940600, 98, 1., 47.13, k_factor=1.),#NLO from MCFM
        sample('ZZ',   'ZZ', '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/AODSIM', 989312, 94, 1.,16.523, k_factor=1.),#NLO from MCFM
       
        sample('WWinclusive', 'WWinclusive','/WWTo2L2Nu_13TeV-powheg/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 1996600, 208, 1., 12.178, k_factor=1.),#already NNLO xs        
        sample('WW200to600', 'WW200to600','/WWTo2L2Nu_Mll_200To600_13TeV-powheg/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 200000, 208, 1., 1.385, k_factor=1.),#already NNLO xs
        sample('WW600to1200', 'WW600to1200','/WWTo2L2Nu_Mll_600To1200_13TeV-powheg/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 200000, 208, 1., 0.0566, k_factor=1.),#already NNLO xs
        sample('WW1200to2500', 'WW1200to2500','/WWTo2L2Nu_Mll_1200To2500_13TeV-powheg/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/AODSIM', 200000, 208, 1., 0.0035, k_factor=1.),#already NNLO xs
        sample('WW2500', 'WW2500','/WWTo2L2Nu_Mll_2500ToInf_13TeV-powheg/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v4/AODSIM', 38969, 208, 1., 0.00005, k_factor=1.),#already NNLO xs
  		
  		sample('dyInclusive50',          'DYInclusive50', '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 28696958, 209 , 1., 6025.2,    k_factor=1., is_madgraph=True),  #(I will update the number of events)

        sample('Wjets', 'Wjets', '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/AODSIM',18995071,52,1.,61526.7,k_factor=1),#already NNLO xs     
        
#        sample('ttbar_pow',     't#bar{t}', '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext4-v1/AODSIM', 186313000, 4 , 1., 815.96, k_factor=1.),
        sample('ttbar_lep',     'ttbar_lep', '/TTTo2L2Nu_13TeV-powheg/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/AODSIM', 104649474, 4 , 1., 87.31, k_factor=1.),

		sample('Wantitop', 'WantiTop', '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',929046,63 , 1., 35.6, k_factor=1.),#already NNLO xs
        sample('tW',     'tW', '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',988400,66 , 1., 35.6, k_factor=1.),#already NNLO xs
        
        
        sample('qcd80to120', 'QCD80to120', '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',6953590,43,1.,2762530,k_factor=1),
        sample('qcd120to170', 'QCD120to170', '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',6848223,43,1.,471100,k_factor=1),
        sample('qcd170to300', 'QCD170to300', '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',6918748,43,1.,117276,k_factor=1),
        sample('qcd300to470', 'QCD300to470', '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',5968960,43,1.,7823,k_factor=1),
        sample('qcd470to600', 'QCD470to600', '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',3977770,43,1.,648.2,k_factor=1),
        sample('qcd600to800', 'QCD600to800', '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',3979884,43,1.,186.9,k_factor=1),
        sample('qcd800to1000', 'QCD800to1000', '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',3973224,43,10,32.293,k_factor=1),
        sample('qcd1000to1400', 'QCD1000to1400', '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',2967947,43,1.,9.4183,k_factor=1),
        sample('qcd1400to1800', 'QCD1400to1800', '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',395725,43,1.,0.84265,k_factor=1),
        sample('qcd1800to2400', 'QCD1800to2400', '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',393760,43,1.,0.114943,k_factor=1),
        sample('qcd2400to3200', 'QCD2400to3200', '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',398452,43,1.,0.00682981,k_factor=1),
        sample('qcd3200', 'QCD3200', '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',391108,43,1.,0.000165445,k_factor=1),

 		#sample('qcd50to80', 'QCD50to80', '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM',50000,43,1.,19204300,k_factor=1),


    ]

samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name
    #if '_c' in sample.name:
    #if 'Zprime' in sample.dataset:
#    sample.ana_dataset = '/%s/federica-%s-f646da20575c2cb2b2eda7b4413fb91e/USER'  % (sample.dataset.split('/')[1], sample.name)
#    if sample.name == 'dy100-200':
#        sample.ana_dataset = '/%s/federica-%s-7a0d7047a2104d11a44d5593620f154b/USER'% (sample.dataset.split('/')[1], sample.name)
#    elif sample.name == 'dy200-400' or sample.name == 'dy400-500' or sample.name == 'dy500-700' or sample.name == 'dy700-800':
#        sample.ana_dataset = '/%s/federica-%s-f646da20575c2cb2b2eda7b4413fb91e/USER'% (sample.dataset.split('/')[1], sample.name)
#    else:

#     
    if 'ZToMuMu' in sample.dataset:
        sample.ana_dataset = '/%s/rradogna-datamc_%s-20330eea4b39a6d27baf680b4cd56b47/USER'  % (sample.dataset.split('/')[1], sample.name)
    if 'W' in sample.name:
     	sample.ana_dataset = '/%s/alfloren-datamc_%s-20330eea4b39a6d27baf680b4cd56b47/USER'  % (sample.dataset.split('/')[1], sample.name)
# 
dyInclusive50.ana_dataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy50-8e76cd6a839f33945c4e1044b64b0194/USER'
ZZ.ana_dataset = '/ZZ_TuneCUETP8M1_13TeV-pythia8/alfloren-datamc_zz-20330eea4b39a6d27baf680b4cd56b47/USER'
Wjets.ana_dataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/alfloren-datamc_wjets-20330eea4b39a6d27baf680b4cd56b47/USER'
WZ.ana_dataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/alfloren-datamc_wz-20330eea4b39a6d27baf680b4cd56b47/USER'
ttbar_lep.ana_dataset = '/TTTo2L2Nu_13TeV-powheg/rradogna-datamc_ttbar_lep-8e76cd6a839f33945c4e1044b64b0194/USER'
# 
WWinclusive.ana_dataset = '/WWTo2L2Nu_13TeV-powheg/alfloren-datamc_WWinclusive-20330eea4b39a6d27baf680b4cd56b47/USER'
# 
qcd80to120.ana_dataset = '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd80to120_PAT_76X_2-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd120to170.ana_dataset = '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd120to170_PAT_76X_2-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd170to300.ana_dataset = '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd170to300_PAT_76X_2-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd300to470.ana_dataset = '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd300to470_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd470to600.ana_dataset = '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd470to600_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd600to800.ana_dataset = '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd600to800_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd800to1000.ana_dataset = '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd800to1000_PAT_76X_2-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd1000to1400.ana_dataset = '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd1000to1400_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd1400to1800.ana_dataset = '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd1400to1800_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd1800to2400.ana_dataset = '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd1800to2400_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd2400to3200.ana_dataset = '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd2400to3200_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'
qcd3200.ana_dataset = '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/cschnaib-datamc_qcd3200_PAT_76X-ae68241ee2858f3bbe9208ae0a97d17f/USER'

#ttbar_pow.ana_dataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/rradogna-datamc_ttbar_pow-8e76cd6a839f33945c4e1044b64b0194/USER'
# #qcd50to80.ana_dataset = ''

# 



__all__ = ['samples'] + [s.name for s in samples]


if __name__ == '__main__':
    if False:
        from dbstools import dbsparents
        for s in samples:
            print s.dataset
            parents = dbsparents(s.dataset)
            for parent in parents:
                for line in os.popen('dbss rel %s' % parent):
                    if 'CMSSW' in line:
                        print parent, line,
            print

    if False:
        import os
        from dbstools import dbsparents
        for s in [ww,wz,zz]:
            print s.dataset
            parents = dbsparents(s.dataset)
            print parents
            os.system('dbsconfig %s > %s' % (parents[-1], s.name))

        os.system('dbss nevents %s' % x.replace('RECO','RAW'))
        os.system('dbss nevents %s' % x)

    if False:
        import os
        from dbstools import dbsparents
        for s in samples:
            print s.dataset
            def fuf(y):
                x = os.popen(y).read()
                for line in x.split('\n'):
                    try:
                        print int(line)
                    except ValueError:
                        pass
            fuf('dbss nevents %s' % s.dataset)
            fuf('dbss nevents %s' % s.dataset.replace('AODSIM','GEN-SIM-RECO'))

    if False:
        for s in samples:
            print s.name
            os.system('grep "total events" ~/nobackup/crab_dirs/384p3/publish_logs/publish.crab_datamc_%s' % s.name)
            os.system('grep "total events" ~/nobackup/crab_dirs/413p2/publish_logs/publish.crab_datamc_%s' % s.name)
            print

    if False:
        os.system('mkdir ~/scratch/wjets')
        for fn in wjets.filenames:
            assert fn.startswith('/store')
            fn = '/pnfs/cms/WAX/11' + fn
            cmd = 'dccp %s ~/scratch/wjets/' % fn
            print cmd
            os.system(cmd)

    if False:
        for s in samples:
            print s.name
            os.system('dbss site %s' % s.dataset)
            print

    if False:
        for s in samples:
            if s.ana_dataset is None:
                continue
            c = []
            for line in os.popen('dbss ana02 find file.numevents where dataset=%s' % s.ana_dataset):
                try:
                    n = int(line)
                except ValueError:
                    continue
                c.append(n)
            c.sort()
            print s.name, c
