from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'StarLightJpsiPhase2_PrivateMC_Step2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step2.py' # The config file for DIGI-L1-DIGI2RAW
config.JobType.maxMemoryMB = 10000             # DIGI step for Hydjet is memory intensive
config.JobType.numCores = 4

config.Data.inputDataset = '/StarLightJpsiPhase2_PrivateMC/phys_heavyions-Step1_CMSSW_14_0_6v2-24797504def601b91dafdc41d7d2f39b/USER'
config.Data.inputDBS = 'phys03'              # Required for USER datasets
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1                  # Process 1 GEN-SIM file per job
config.Data.publication = True
config.Data.outputDatasetTag = 'Step2_DIGI_RAW_CMSSW_14_0_6v2'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/davidlw/' 

config.Site.storageSite = 'T2_CH_CERN'
