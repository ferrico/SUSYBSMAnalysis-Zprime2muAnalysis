import FWCore.ParameterSet.Config as cms

ResolutionAtZ = cms.EDAnalyzer('ResolutionAtZ',
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   )

ResolutionAtZ_MiniAOD = cms.EDAnalyzer('ResolutionAtZ',
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   )
