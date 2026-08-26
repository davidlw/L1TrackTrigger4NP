# ---------------------------------------------------------------------------
# crab_offlineTrackNtuple_cfg.py
#
# Run offlineTrackNtuple_cfg.py over the published step3 datasets.
#
#   1. pick a sample below
#   2. cmsenv in the release where OfflineTrackNtupleMaker is BUILT
#   3. crab submit -c crab_offlineTrackNtuple_cfg.py
#
# IMPORTANT -- submit from CMSSW_15_1_0_patch3, not CMSSW_14_0_6.
# CRAB ships the local release area (lib/, python/) with the job, so the plugin
# has to exist in the area you submit from. The samples were *produced* in
# CMSSW_14_0_6, but reading them from 15_1_0_patch3 is fine and is what was
# tested. Copy this file and offlineTrackNtuple_cfg.py next to each other in
# CMSSW_15_1_0_patch3/src and submit from there.
#
# The dataset names below are the *output* datasets of the Step3 crab tasks in
# CMSSW_14_0_6/src/crab_projects (the tag reads "Step2_RECO-MINIAOD..." but these
# are the step-3 RECO-MINIAOD outputs). All five were produced with the extra
# keeps, so the offline products the matching needs are present.
# ---------------------------------------------------------------------------

from CRABClient.UserUtilities import config
config = config()

# --- pick one -------------------------------------------------------------
SAMPLE = 'QEDMuMu'

# unitsPerJob is files-per-job. The UPC samples are small and quick; HYDJET PbPb
# and EPOS pPb have far higher multiplicity, so the hit-based association costs
# much more time and memory -- one file per job and a bigger memory request.
SAMPLES = {
    'QEDMuMu': dict(
        dataset='/StarLightQEDMuMuPhase2_PrivateMC/phys_heavyions-Step2_RECO-MINIAOD_CMSSW_14_0_6_v4-51f0e818e6fa09b9e74bf4ad8ded9703/USER',
        unitsPerJob=5, memoryMB=2500),
    'QEDEE': dict(
        dataset='/StarLightQEDEEPhase2_PrivateMC/phys_heavyions-Step2_RECO-MINIAOD_CMSSW_14_0_6_v4-51f0e818e6fa09b9e74bf4ad8ded9703/USER',
        unitsPerJob=5, memoryMB=2500),
    'JpsiMuMu': dict(
        dataset='/StarLightJpsiMuMuPhase2_PrivateMC/phys_heavyions-Step2_RECO-MINIAOD_CMSSW_14_0_6_v4-51f0e818e6fa09b9e74bf4ad8ded9703/USER',
        unitsPerJob=5, memoryMB=2500),
    'HydjetPbPb': dict(
        dataset='/HydjetPbPbPhase2_PrivateMC/phys_heavyions-Step2_RECO-MINIAOD_CMSSW_14_0_6_v4-51f0e818e6fa09b9e74bf4ad8ded9703/USER',
        unitsPerJob=1, memoryMB=8000),
    'EPOSpPb': dict(
        dataset='/EPOSpPbPhase2_PrivateMC/phys_heavyions-Step2_RECO-MINIAOD_CMSSW_14_0_6_v4-51f0e818e6fa09b9e74bf4ad8ded9703/USER',
        unitsPerJob=1, memoryMB=6000),
}
S = SAMPLES[SAMPLE]

config.General.requestName = 'OfflineTrackNtuple_%s_v1' % SAMPLE
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'offlineTrackNtuple_cfg.py'
config.JobType.maxMemoryMB = S['memoryMB']
# Must equal process.options.numberOfThreads in the pset. offlineTrackNtuple_cfg.py
# does not set it, so it is 1; CRAB refuses the task if these disagree.
config.JobType.numCores = 1

config.Data.inputDataset = S['dataset']
config.Data.inputDBS = 'phys03'          # required for USER datasets
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = S['unitsPerJob']
config.Data.publication = False          # plain ntuples, nothing to publish
config.Data.outputDatasetTag = 'OfflineTrackNtuple_%s_v1' % SAMPLE
config.Data.outLFNDirBase = '/store/group/phys_heavyions/davidlw/'

config.Site.storageSite = 'T2_CH_CERN'
