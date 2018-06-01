
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_effres_dy4500to6000'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'
#config.JobType.priority = 1

config.Data.inputDataset =  '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 20000
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_dy4500to6000'
config.Data.outLFNDirBase = '/store/user/ferrico'
config.Data.ignoreLocality = False
# config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_IT_Bari'
                          
