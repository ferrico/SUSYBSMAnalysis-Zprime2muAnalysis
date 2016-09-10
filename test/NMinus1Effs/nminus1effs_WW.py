#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff import loose_cut, trigger_match, tight_cut, allDimuons

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels, rec_level_module
#tracks = ['global', 'inner', 'tpfms', 'picky', 'tunep', 'tmr', 'tunepnew']
tracks = ['tunepnew', 'inner']
rec_levels(process, tracks)

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source.fileNames =[#'file:PAT_SingleMuRun2015B-Rereco-Suite_251162_251559_20160120153115/crab_SingleMuRun2015B-Rereco-Suite_251162_251559_20160120153115/results/Zprime_123.root',
                          'file:root://xrootd-cms.infn.it///store/user/alfloren/PAATuples/WWTo2L2Nu_13TeV-powheg/datamc_WWinclusive/160524_124813/0000/pat_1.root'                          
#'file:root://xrootd-cms.infn.it///store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/datamc_dy50to120/160509_211446/0000/pat_1.root',
#                         'file:root://xrootd-cms.infn.it///store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/datamc_dy50to120/160509_211446/0000/pat_2.root', 
                        #'file:root://xrootd-cms.infn.it///store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/datamc_dy50to120/160509_211446/0000/pat_3.root',
#       'file:root://xrootd-cms.infn.it///store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/datamc_dy50to120/160509_211446/0000/pat_4.root',                         
#root://xrootd-cms.infn.it///store/mc/RunIISpring16DR80/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/04A54643-BA02-E611-8778-00266CFCC618.root',
                           ]
#process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
#readFiles.extend( [
#                          'file:root://xrootd-cms.infn.it///store/user/rradogna/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/datamc_dy50to120/160509_211446/0000/pat_1.root',
#'root://xrootd-cms.infn.it///store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/datamc_dy120to200/160509_211431/0000/pat_1.root'
#/store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/datamc_dy1400to2300/160509_211325/0000/pat_1.root',
#                   '/store/user/rradogna/RelValZMM_13/datamc_dy50_startup/150730_143145/0000/pat_1.root'
#                   ] );


#secFiles.extend( [
#               ] )

process.maxEvents.input = -1

# Define the numerators and denominators, removing cuts from the
# allDimuons maker. "NoX" means remove cut X entirely (i.e. the
# loose_cut denominators), "TiX" means move cut X from the loose_cut
# to the tight_cut (meaning only one muon instead of two has to pass
# the cut).  "NoNo" means remove nothing (i.e. the numerator). This
# will break if loose_, tight_cut strings are changed upstream, so we
# try to check those with a simple string test below.

process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
                                       src = cms.InputTag('prunedMCLeptons'),
                                       min_mass = cms.double(50),
                                       max_mass = cms.double(200),
                                       )


cuts = [
    ('Pt',      'pt > 53'),
    ('DB',      'abs(dB) < 0.2'),
    ('Iso',     'isolationR03.sumPt / innerTrack.pt < 0.10'),
    ('TkLayers','globalTrack.hitPattern.trackerLayersWithMeasurement > 5'),
    ('PxHits',  'globalTrack.hitPattern.numberOfValidPixelHits >= 1'),
    ('MuHits',  'globalTrack.hitPattern.numberOfValidMuonHits > 0'),
    ('MuMatch', ('numberOfMatchedStations > 1', 'isTrackerMuon')),
    ]

for name, cut in cuts:
    if type(cut) != tuple:
        cut = (cut,)

    lc = loose_cut
    for c in cut:
        if c not in lc:
            raise ValueError('cut "%s" not in cut string "%s"' % (c, lc))
        lc = lc.replace(' && ' + c, '') # Relies on none of the cuts above being first in the list.

    obj_no = allDimuons.clone(loose_cut = lc)
#    obj_no = allDimuons.clone(loose_cut = lc,tight_cut = tight_cut.replace(trigger_match, ''))#N-2
    setattr(process, 'allDimuonsNo' + name, obj_no)
    
    obj_ti = obj_no.clone(tight_cut = tight_cut + ' && ' + ' && '.join(cut))
    setattr(process, 'allDimuonsTi' + name, obj_ti)

process.allDimuonsNoNo      = allDimuons.clone()
#process.allDimuonsNoNo      = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))#N-2
process.allDimuonsNoTrgMtch = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))

alldimus = [x for x in dir(process) if 'allDimuonsNo' in x or 'allDimuonsTi' in x]

# Sanity check that the replaces above did something.
for x in alldimus:
    if 'NoNo' in x:
        continue
    o = getattr(process, x)
    assert o.loose_cut.value() != loose_cut or o.tight_cut.value() != tight_cut

process.p = cms.Path(process.goodDataFilter * process.DYGenMassFilter * process.muonPhotonMatch * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))

# For all the allDimuons producers, make dimuons producers, and
# analyzers to make the histograms.
for alld in alldimus:
    dimu = process.dimuons.clone(src = alld)
    name = alld.replace('allD', 'd')
    setattr(process, name, dimu)
    hists = HistosFromPAT.clone(dilepton_src = name, leptonsFromDileptons = True)
    setattr(process, name.replace('dimuons', ''), hists)
    process.p *= dimu * hists

# Handle the cuts that have to be applied at the
# Zprime2muCompositeCandidatePicker level.
#process.allDimuonsN2 = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))#N-2
#process.p *= process.allDimuonsN2#N-2
process.dimuonsNoB2B     = process.dimuons.clone()
process.dimuonsNoVtxProb = process.dimuons.clone()
process.dimuonsNoDptPt   = process.dimuons.clone()
#process.dimuonsNoB2B     = process.dimuons.clone(src = 'allDimuonsN2')#N-2
#process.dimuonsNoVtxProb = process.dimuons.clone(src = 'allDimuonsN2')#N-2
#process.dimuonsNoDptPt   = process.dimuons.clone(src = 'allDimuonsN2')#N-2
delattr(process.dimuonsNoB2B,     'back_to_back_cos_angle_min')
delattr(process.dimuonsNoVtxProb, 'vertex_chi2_max')
delattr(process.dimuonsNoDptPt,   'dpt_over_pt_max')
process.p *= process.allDimuons
for dimu in ['dimuonsNoB2B', 'dimuonsNoVtxProb', 'dimuonsNoDptPt']:
    hists = HistosFromPAT.clone(dilepton_src = dimu, leptonsFromDileptons = True)
    setattr(process, dimu.replace('dimuons', ''), hists)
    process.p *= getattr(process, dimu) * hists

# Special case to remove |dB| and B2B cuts simultaneously, as they can
# be correlated (anti-cosmics).
process.allDimuonsNoCosm = process.allDimuons.clone(loose_cut = loose_cut.replace(' && abs(dB) < 0.2', ''))
#process.allDimuonsNoCosm = process.allDimuons.clone(loose_cut = loose_cut.replace(' && abs(dB) < 0.2', ''), tight_cut = tight_cut.replace(trigger_match, '')) #N-2
process.dimuonsNoCosm = process.dimuons.clone(src = 'allDimuonsNoCosm')
delattr(process.dimuonsNoCosm, 'back_to_back_cos_angle_min')
process.NoCosm = HistosFromPAT.clone(dilepton_src = 'dimuonsNoCosm', leptonsFromDileptons = True)
process.p *= process.allDimuonsNoCosm * process.dimuonsNoCosm * process.NoCosm


process.p *= rec_level_module(process, HistosFromPAT,     'Histos',     tracks)
#process.p *= rec_level_module(process, process.ResolutionUsingMC, 'Resolution', tracks)

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_nminus1_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nminus1effs_WW.py'
#config.JobType.priority = 1

config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'phys03'
job_control
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/ferrico'

config.Site.storageSite = 'T2_IT_Bari'
#config.Site.storageSite = 'T2_IT_Legnaro'

'''

    just_testing = 'testing' in sys.argv
    if not 'no_data' in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Run2016MuonsOnly_ll
        Run2016MuonsOnly_ll.writeJSON('tmp.json')

        dataset_details = [
#            ('SingleMuonRun2015B-Prompt_251162_251499',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251162_251499_20150713100409-3aa7688518cb1f1b044caf15b1a9ed05/USER'),
#            ('SingleMuonRun2015B-Prompt_251500_251603',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251500_251603_20150718235715-9996471c14459acaec01707975d1e954/USER'),
#            ('SingleMuonRun2015B-Prompt_251613_251883',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251613_251883_20150719000207-9996471c14459acaec01707975d1e954/USER'),                           
#            ('SingleMuonRun2015C-Prompt_253888_254914',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015C-Prompt_253888_254914_20150831150018-681693e882ba0f43234b3b41b1bbc39d/USER'),
#			 ('SingleMuonRun2015C25ns-rereco76_254227_254907', '/SingleMuon/alfloren-SingleMuRun2015C25ns-Rereco_254227_254907_20160120153639-332cf72ab044858cbe7c1d1b03f22dbc/USER'),
#             ('SingleMuonRun2015D-rereco76_256630_260627_cschnaib', '/SingleMuon/cschnaib-datamc_SingleMuRun2015D-Rereco76X_PAT-843ac0dcce157982e3f7d22621d7dc4b/USER'),
            ('SingleMuonRun2016B-Prompt-v2_273150_273730_rradogna', '/SingleMuon/rradogna-datamc_SingleMuonRun2016B-Prompt-v2_273150_273730_20160530153025-02d6fdda0cfc6c7b5229d43ba172d9c1/USER'),
            ('SingleMuonRun2016B-Prompt-v2_273731_274421_rradogna', '/SingleMuon/rradogna-datamc_SingleMuonRun2016B-Prompt-v2_273731_274421_20160612232019-02d6fdda0cfc6c7b5229d43ba172d9c1/USER'),

            ]

        for name, ana_dataset in dataset_details:
            print name

            new_py = open('nminus1effs.py').read()
            new_py += "\nprocess.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'\n"
            open('nminus1effs_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()
            job_control = '''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'tmp.json' #######
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_254833_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274421_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
'''
            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                #os.system('crab submit -c crabConfig.py --dryrun') #--dryrun
                os.system('crab submit -c crabConfig.py') #--dryrun

        if not just_testing:
            os.system('rm crabConfig.py nminus1effs_crab.py nminus1effs_crab.pyc tmp.json')

    if not 'no_mc' in sys.argv:
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 10000
''')

        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
        #samples =[dy50to120,DY120to200Powheg,DY200to400Powheg,DY400to800Powheg,DY800to1400Powheg,dy1400to2300, ttbar_pow]#,ttbar, wz, ww_incl, zz_incl, dy50to120]
        #samples =[dy50, dy120, dy200, dy400, dy800, dy1400, dy2300, dy3500, dy4500, dy6000, dy7500, dy8500, dy9500, zpsi5000, ttbar, inclmu15]
        
        
        ########################

        
        #tutti i samples
        samples =[WWinclusive]
		#dyInclusive50,WW200to600, WW600to1200, WW1200to2500, WW2500,

		########################


		#solo QCD
        #samples =[qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200]


		########################

        for sample in samples:
            print sample.name
            open('crabConfig.py', 'wt').write(crab_cfg % sample)
            if not just_testing:
                #os.system('crab submit -c crabConfig.py --dryrun')
				os.system('crab submit -c crabConfig.py')
    
        if not just_testing:
            os.system('rm crabConfig.py')
