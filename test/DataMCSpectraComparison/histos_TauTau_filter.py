#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
#from SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi import DYGenMassFilter



process.source.fileNames =[
                          'file:root://xrootd-cms.infn.it///store/user/rradogna/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/datamc_dy50/160613_110131/0000/pat_75.root',
                           ]
process.maxEvents.input =50
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'
#process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True

# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
#import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionOld_cff as OurSelectionOld
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2011EPS_cff as OurSelection2011EPS
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff as OurSelectionNew
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff as OurSelectionDec2012

# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
dils = [
    ('MuonsPlusMuonsMinusBarrel',    '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0 && daughter(0).eta()<=1.2 && daughter(1).eta()<=1.2 && daughter(0).eta()>=-1.2 &&  daughter(1).eta()>=-1.2'),
    ('MuonsPlusMuonsMinusPositive',    '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0 && ((daughter(0).eta()>1.2 && daughter(1).eta()>-1.2) || (daughter(1).eta()>1.2 && daughter(0).eta()>-1.2))'),
    ('MuonsPlusMuonsMinusNegative',    '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0 && (daughter(0).eta()<=-1.2 || daughter(1).eta()<=-1.2)'),
    ('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    ('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
    ('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
    ('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
    ('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
    ('MuonsPlusElectronsMinus',      '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == -2'),
    ('MuonsMinusElectronsPlus',      '%(leptons_name)s:muons@- %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == 2'),
    ('MuonsPlusElectronsPlus',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == -24'),
    ('MuonsMinusElectronsMinus',     '%(leptons_name)s:muons@- %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == 24'),
    ('MuonsElectronsOppSign',        '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     ''),
    ('MuonsElectronsSameSign',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
    ('MuonsElectronsAllSigns',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
    ]

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
#    'VBTF'     : VBTFSelection,
#    'OurOld'   : OurSelectionOld,
#    'OurEPS'   : OurSelection2011EPS,
    #'OurNew'   : OurSelectionNew,

   # 'EtaBarrel' : OurSelectionDec2012,
#    'EtaEpositive' : OurSelectionDec2012,
#    'EtaNegative' : OurSelectionDec2012,
    'Our2012'  : OurSelectionDec2012,
    #'OurNoIso' : OurSelectionDec2012,
    #'EmuVeto'  : OurSelectionDec2012,
    'Simple'   : OurSelectionDec2012, # The selection cuts in the module listed here are ignored below.
#    'VBTFMuPrescaled' : VBTFSelection,
    #'OurMuPrescaledNew'  : OurSelectionNew,
  # 'OurMuPrescaled2012' : OurSelectionDec2012
    }


process.DYGenMassFilter = cms.EDFilter('TauTauSelection',
                                       src = cms.InputTag('prunedMCLeptons'),
                                       min_mass = cms.double(50),
                                       max_mass = cms.double(200),
                                       )

# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for cut_name, Selection in cuts.iteritems():
    # Keep track of modules to put in the path for this set of cuts.
    path_list = []
   
    
    #path_list.append(DYGenMassFilter)

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
    # just marks each muon with a userInt "cutFor" that is 0 if it
    # passes the cuts, and non-0 otherwise; it does not actually drop
    # any of the muons. The cutFor flag actually gets ignored by the
    # LooseTightPairSelector in use for all the cuts above, at
    # present.
    leptons_name = cut_name + 'Leptons'
    DYGenMassFilter_name = cut_name + 'GenFilter'

    if cut_name == 'Simple':
        muon_cuts = ''
   # elif 'MuPrescaled' in cut_name:
      #  muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut
    leptons = process.leptons.clone(muon_cuts = muon_cuts)
    DYGenMassFilter = process.DYGenMassFilter.clone()
    
    if cut_name == 'EmuVeto':
        leptons.electron_muon_veto_dR = 0.1
    # Keep using old TuneP for past selections
    if 'Dec2012' not in Selection.__file__:
        leptons.muon_track_for_momentum = cms.string('TuneP')
    setattr(process, leptons_name, leptons)
    setattr(process, DYGenMassFilter_name, DYGenMassFilter)
    
    path_list.append(DYGenMassFilter*leptons)
     #path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        # For the EmuVeto path, we only care about e-mu events.
        if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
            continue

        # For the MuPrescaled paths, we don't care about e-mu events.
        if 'MuPrescaled' in cut_name and 'Electron' in dil_name:
            continue

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)
        dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        # Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        if cut_name == 'Simple':
            alldil.electron_cut_mask = cms.uint32(0)
            #alldil.loose_cut = 'isGlobalMuon && pt > 20.'#to be changed for first runs
            alldil.loose_cut = 'isGlobalMuon && pt > 20.'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.sort_by_pt = True
            dil.do_remove_overlap = False
            if hasattr(dil, 'back_to_back_cos_angle_min'):
                delattr(dil, 'back_to_back_cos_angle_min')
            if hasattr(dil, 'vertex_chi2_max'):
                delattr(dil, 'vertex_chi2_max')
            if hasattr(dil, 'dpt_over_pt_max'):
                delattr(dil, 'dpt_over_pt_max')
        elif cut_name == 'OurNoIso':
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', '')
       # elif 'MuPrescaled' in cut_name:
          #  alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
         #   assert alldil.tight_cut == trigger_match
         #   alldil.tight_cut = prescaled_trigger_match

        # Histos now just needs to know which leptons and dileptons to use.
        histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
       
        path_list.append(alldil * dil * histos)

    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    pobj = process.muonPhotonMatch * reduce(lambda x,y: x*y, path_list)
    if 'VBTF' not in cut_name and cut_name != 'Simple':
        pobj = process.goodDataFilter * pobj
    if 'MuPrescaled' in cut_name: ####### Now it seams that there are no prescaled path ########
        pobj = process.PrescaleToCommon * pobj ####### Now it seams that there are no prescaled path ########
    path = cms.Path(pobj)
    setattr(process, pathname, path)
   


def ntuplify(process, fill_gen_info=False):
    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler',
                                           dimu_src = cms.InputTag('SimpleMuonsAllSigns'),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlinePrimaryVertices'),
					   TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),
                                           genEventInfo = cms.untracked.InputTag('generator')
                                           )
    process.SimpleNtuplerEmu = process.SimpleNtupler.clone(dimu_src = cms.InputTag('SimpleMuonsElectronsAllSigns'))

    if fill_gen_info:
        from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
        process.SimpleNtupler.hardInteraction = hardInteraction
        
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.SimpleNtupler * process.SimpleNtuplerEmu
ntuplify(process) #to have ntuples also running in interactive way

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')

    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD','','HLT')
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.triggerSummaryAnalyzerAOD

    process.PrintOriginalMuons = cms.EDAnalyzer('PrintEvent', muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
    process.pathSimple *= process.PrintOriginalMuons

    pe = process.PrintEventSimple = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.PrintEventSimple

    #- 2011-2012 selection (Nlayers > 8)
    #process.PrintEventOurNew = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
    #process.PrintEventOurNewSS = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
    #process.PrintEventOurNewEmu = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
    #process.pathOurNew *= process.PrintEventOurNew * process.PrintEventOurNewSS * process.PrintEventOurNewEmu

    #- December 2012 selection (Nlayers > 5, re-tuned TuneP, dpT/pT < 0.3)
    if hasattr(process, 'pathOur2012'):
        process.PrintEventOur2012    = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsPlusMuonsMinus'))
        process.PrintEventOur2012SS  = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsSameSign'))
        process.PrintEventOur2012Emu = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsElectronsOppSign'))
        process.pathOur2012 *= process.PrintEventOur2012 * process.PrintEventOur2012SS * process.PrintEventOur2012Emu

def check_prescale(process, trigger_paths, hlt_process_name='HLT'):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
    process.CheckPrescale.trigger_paths = cms.vstring(*trigger_paths)
    process.pCheckPrescale = cms.Path(process.CheckPrescale)

def for_data(process):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'	
    ntuplify(process)
    #check_prescale(process, trigger_paths) ####### Now it seams that there are no prescaled path ########

def for_mc(process, hlt_process_name, fill_gen_info):
    ntuplify(process, fill_gen_info)
    switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)

def get_dataset(run):
    #JMTBAD common with dataset_details in submit below, make a DataSamples.py?
    run = int(run)
    if 190450 <= run <= 191284:
        return '/SingleMu/tucker-datamc_SingleMuRun2012A_Prompt_190450_191284_20120418134612-57b19813ab8f2ab142c4566dc6738156/USER'
    else:
        raise ValueError('dunno how to do run %i' % run)

if 'int_data' in sys.argv:
    for_data(process)
    printify(process)
    
if 'int_mc' in sys.argv:
    for_mc(process, 'HLT', False)
    printify(process)
    
if 'gogo' in sys.argv:
    for_data(process)
    printify(process)
    
    n = sys.argv.index('gogo')
    run, lumi, event = sys.argv[n+1], sys.argv[n+2], sys.argv[n+3]
    print run, lumi, event
    run = int(run)
    lumi = int(lumi)
    event = int(event)
    filename = [x for x in sys.argv if x.endswith('.root')]
    if filename:
        filename = filename[0]
    else:
        dataset = get_dataset(run)
        print dataset
        output = os.popen('dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query="find file where dataset=%s and run=%s and lumi=%s"' % (dataset, run, lumi)).read()
        print repr(output)
        filename = [x for x in output.split('\n') if x.endswith('.root')][0]
    print filename
    process.source.fileNames = [filename]
    from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process
    set_events_to_process(process, [(run, event)])

f = file('outfile', 'w')
f.write(process.dumpPython())
f.close()

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_datamc_%(name)s_JsonPromptReco'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histosWW.py'   
#config.JobType.priority = 1

config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'phys03'
job_control
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/ferrico'
config.Data.ignoreLocality = True 

config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_IT_Bari'


'''
    
    just_testing = 'testing' in sys.argv
        
    # Run on data.
    if 'no_data' not in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *

        dataset_details = [
            # ('SingleMuRun2015D-Express_256843_257490',    '/ExpressPhysics/alfloren-SingleMuRun2015D-Express_256584_257490_20150927202323-fb80d99301269fbfefa76b46bb4235ae/USER'),
             # ('SingleMuonRun2015D-Prompt_256629_256842',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015D-Prompt_256629_256842_20150926113604-c9b39dd88dc98b683a1d7cecc8f6c42c/USER'),
             # ('SingleMuonRun2015D-Prompt_256843_257819',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015D-Prompt_256843_257819_20151002140028-c9b39dd88dc98b683a1d7cecc8f6c42c/USER'),
             # ('SingleMuonRun2015D-Prompt_257820_258157',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015D-Prompt_257820_258157_20151004232715-c9b39dd88dc98b683a1d7cecc8f6c42c/USER')
                           #('SingleMuonRun2015D-Prompt_258158_258158',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015D-Prompt_258158_258158_20151009194535-c9b39dd88dc98b683a1d7cecc8f6c42c/USER'),
                           #('SingleMuonRun2015D-Prompt_258159_258432',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015D-Prompt_258159_258432_20151009203706-c9b39dd88dc98b683a1d7cecc8f6c42c/USER'),
                # ('SingleMuonRun2015B-rereco76_251028_251160',    '/SingleMu/alfloren-SingleMuRun2015B-Rereco_251028_251160_20160120151841-332cf72ab044858cbe7c1d1b03f22dbc/USER'),
                #('SingleMuonRun2015B-rereco76_251162_251559',    '/SingleMuon/alfloren-SingleMuRun2015B-Rereco-Suite_251162_251559_20160120153115-332cf72ab044858cbe7c1d1b03f22dbc/USER'),
                #('SingleMuonRun2015C50ns-rereco76_254883_255899', '/SingleMuon/alfloren-SingleMuRun2015C50ns-Rereco_254883_255899_20160120153444-332cf72ab044858cbe7c1d1b03f22dbc/USER'),
                 ('SingleMuonRun2015C25ns-rereco76_254227_254907', '/SingleMuon/alfloren-SingleMuRun2015C25ns-Rereco_254227_254907_20160120153639-332cf72ab044858cbe7c1d1b03f22dbc/USER'),
                 ('SingleMuonRun2015D-rereco76_256630_260627', '/SingleMuon/rradogna-datamc_SingleMuonRun2015D-Rereco_256630_260627_20160204164737-843ac0dcce157982e3f7d22621d7dc4b/USER'),

            ]

        lumi_lists = [
           # 'NoLumiMask'
  #           'DCSOnly',
#            'Run2012PlusDCSOnlyMuonsOnly',
            'Run2016MuonsOnly',
           # 'Run2015',
            ]

        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll') if lumi_name != 'NoLumiMask' else None
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))
                
        for dataset_name, ana_dataset, lumi_name, lumi_list in jobs:
            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = json_fn

            name = '%s_%s' % (lumi_name, dataset_name)
            print name

            new_py = open('histos.py').read()
            new_py += "\nfor_data(process)\n"
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
config.Data.splitting = 'LumiBased'
#config.Data.runRange = '256843-257490'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 200
config.Data.lumiMask = '%(lumi_mask)s' #######
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 10000
    ''')

        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples

        combine_dy_samples = len([x for x in samples if x.name in ['dy200', 'dy400', 'dy500', 'dy700', 'dy800', 'dy1500', 'dy2000', 'dy3000']]) > 0
        print 'combine_dy_samples:', combine_dy_samples

        for sample in reversed(samples):
            print sample.name

            new_py = open('histos.py').read()
            sample.fill_gen_info = sample.name in ['dy50', 'dy100to200', 'dy200to400', 'dy400', 'dy800', 'dy1400', 'dy2300', 'dy3500', 'dy4500', 'dy6000', 'dy7500', 'dy8500', 'dy9500', 'zpsi5000']
            new_py += "\nfor_mc(process, hlt_process_name='%(hlt_process_name)s', fill_gen_info=%(fill_gen_info)s)\n" % sample

            if combine_dy_samples and (sample.name == 'zmumu' or 'dy' in sample.name):
                mass_limits = {
                    #'dy50'      : (  50,     120),
                    #'dy120'     : ( 120,     200),
                    'dy200'     : ( 100,     200),
                    'dy400'     : ( 200,     400),
                    'dy500'     : ( 400,     500),
                    'dy700'     : ( 500,     700),
                    'dy800'     : ( 700,     800),
                    'dy1500'    : (1000,    1500),
                    'dy2000'    : (1500,    2000),
                    'dy3000'    : (2000,    3000),
                    #'dy7500'    : (6000,    7500),
                    #'dy8500'    : (8500,    9500),
                    #'dy9500'    : (9500,  100000),
                    }
                lo,hi = mass_limits[sample.name]
                from SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi import dy_gen_mass_cut
                new_cut = dy_gen_mass_cut % locals()

                new_py += '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi')

process.DYGenMassFilter.cut = "%(new_cut)s"
for pn,p in process.paths.items():
    setattr(process, pn, cms.Path(process.DYGenMassFilter*p._seq))
''' % locals()

            open('histos_crab.py', 'wt').write(new_py)

            open('crabConfig.py', 'wt').write(crab_cfg % sample)
            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crabConfig.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')