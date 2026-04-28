from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'StarLightJpsiPhase2_PrivateMC_Step1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'step1.py' # The file created by cmsDriver
# Increase memory/time for chained steps
config.JobType.maxMemoryMB = 3000 
config.JobType.numCores = 1

config.Data.outputPrimaryDataset = 'StarLightJpsiPhase2_PrivateMC'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 2000  # Number of events per job
NJOBS = 100  # Total number of jobs
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
config.Data.outputDatasetTag = 'Step1_CMSSW_14_0_6v2'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/davidlw/' 

# Edit the site to your assigned storage site (e.g., T2_US_Purdue)
config.Site.storageSite = 'T2_CH_CERN'
