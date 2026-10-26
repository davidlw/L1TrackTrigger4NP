from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'StarLightQEDEEPhase2_PrivateMC_L1TrackNtuplev2'
#config.General.requestName = 'StarLightQEDMuMuPhase2_PrivateMC_L1TrackNtuplev2'
#config.General.requestName = 'EPOSpPbPhase2_PrivateMC'
config.General.requestName = 'HydjetPbPbPhase2_PrivateMC'
#config.General.requestName = 'MinBias_TuneCP5_14TeV-pythia8'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'rerunL1_trackhitntuple_cfg.py' # The config file for DIGI-L1-DIGI2RAW
config.JobType.maxMemoryMB = 5000             # DIGI step for Hydjet is memory intensive
config.JobType.numCores = 2

#config.Data.inputDataset = '/StarLightQEDEEPhase2_PrivateMC/phys_heavyions-Step2_DIGI_RAW_CMSSW_14_0_6_tpminpt0p3-706485a162e6a5bb49a23e5c839ad1cc/USER' 
#config.Data.inputDataset = '/StarLightQEDMuMuPhase2_PrivateMC/phys_heavyions-Step2_DIGI_RAW_CMSSW_14_0_6_v4_tpminpt0p3-02742f24b5cd022229ea44e62d5e3f48/USER'
#config.Data.inputDataset = '/EPOSpPbPhase2_PrivateMC/phys_heavyions-Step2_DIGI_RAW_CMSSW_14_0_6_tpminpt0p3-0cf1277c45fb77e5f8de2cea235e1e2d/USER' 
config.Data.inputDataset = '/HydjetPbPbPhase2_PrivateMC/phys_heavyions-Step2_DIGI_RAW_CMSSW_14_0_6_tpminpt0p3-6c26596e10aa1561d97e290183c86cef/USER'
#config.Data.inputDataset = '/MinBias_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU140_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD'
config.Data.inputDBS = 'phys03'              # Required for USER datasets
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.unitsPerJob = 5
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_heavyions/nalewis/' 

config.Site.storageSite = 'T2_CH_CERN'
