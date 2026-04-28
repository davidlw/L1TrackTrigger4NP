from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'StarLightJpsiPhase2_PrivateMC_L1TrackNtuplev2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'rerunL1_trackhitntuple_cfg.py' # The config file for DIGI-L1-DIGI2RAW
config.JobType.maxMemoryMB = 5000             # DIGI step for Hydjet is memory intensive
config.JobType.numCores = 2

config.Data.inputDataset = '/StarLightJpsiPhase2_PrivateMC/phys_heavyions-Step2_DIGI_RAW_CMSSW_14_0_6v2-ed967c326098feeffad433ab0d0e9083/USER'
config.Data.inputDBS = 'phys03'              # Required for USER datasets
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10                  # Process 1 GEN-SIM file per job
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_heavyions/davidlw/' 

config.Site.storageSite = 'T2_CH_CERN'
