
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: -s L1TrackTrigger,L1,L1P2GT,NANO:@Phase2L1DPGwithGen --conditions auto:phase2_realistic_T33 --geometry ExtendedRun4D110 --era Phase2C17I13M9 --eventcontent NANOAOD --datatier GEN-SIM-DIGI-RAW-MINIAOD --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110 --filein root://cmsxrootd.fnal.gov///store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000/11d1f6f0-5f03-421e-90c7-b5815197fc85.root --fileout file:output_Phase2_L1T.root --python_filename rerunL1_cfg.py --inputCommands=keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT --mc -n 100 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('NANO',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.SimPhase2L1GlobalTriggerEmulator_cff')
process.load('L1Trigger.Configuration.Phase2GTMenus.SeedDefinitions.step1_2024.l1tGTMenu_cff')
process.load('DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
#    fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov///store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000/11d1f6f0-5f03-421e-90c7-b5815197fc85.root'),
#    fileNames = cms.untracked.vstring('file:step3.root'),
 fileNames = cms.untracked.vstring(
'file:/eos/cms/store/group/phys_heavyions/davidlw/StarLightJpsiPhase2_PrivateMC/Step2_DIGI_RAW_CMSSW_14_0_6/260325_123236/0000/step2_1.root'
#'file:/eos/cms/store/group/phys_heavyions/davidlw/L1TrackTrigger/hydjet/step2_hydjet_1.root',
#'file:/eos/cms/store/group/phys_heavyions/davidlw/L1TrackTrigger/hydjet/step2_hydjet_2.root',
#'file:/eos/cms/store/group/phys_heavyions/davidlw/L1TrackTrigger/hydjet/step2_hydjet_3.root',
#'file:/eos/cms/store/group/phys_heavyions/davidlw/L1TrackTrigger/hydjet/step2_hydjet_4.root',
),
# fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_heavyions/davidlw/L1TrackTrigger/hydjet/step3_hydjet.root'),
# fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov///store/relval/CMSSW_14_1_0_pre3/RelValHydjetQMinBias_5362GeV/GEN-SIM-RECO/140X_mcRun4_realistic_v3_STD_2026D98_HIN_noPU-v1/2590000/9090e79c-cdcd-4606-9362-96f064976a7c.root'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop l1tPFJets_*_*_*',
        'drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(2),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('-s nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Path and EndPath definitions
#process.StubRECO_step = cms.Path(process.TrackTriggerStubs)
#process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000 

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customisePhase2TTOn110
from L1Trigger.Configuration.customisePhase2TTOn110 import customisePhase2TTOn110 

#call to customisation function customisePhase2TTOn110 imported from L1Trigger.Configuration.customisePhase2TTOn110
process = customisePhase2TTOn110(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

'''
# --- FORCE UPC 0.8 GeV RESET ---
def apply_upc_low_pt_settings(process):
# --- MANUAL WINDOW SCALING FOR 0.8 GeV ---
    if hasattr(process, 'TTStubAlgorithm_official_Phase2TrackerDigi_'):
        # 1. Scale the Barrel Cuts
        original_barrel = process.TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut
        process.TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble([x * 7 for x in original_barrel])

        # 2. Scale the Tilted Barrel Cuts
        for pset in process.TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet:
            original_tilted = pset.TiltedCut
            pset.TiltedCut = cms.vdouble([x * 7 for x in original_tilted])

        # 3. Scale the Endcap Cuts
        for pset in process.TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet:
            original_endcap = pset.EndcapCut
            pset.EndcapCut = cms.vdouble([x * 7 for x in original_endcap])

        # Define the Analytic calculation for 0.8 GeV
#        process.TTStubAlgorithm_official_Phase2TrackerDigi_.findingAlgorithm = cms.PSet(
#            AlgorithmName = cms.string("TTStubAlgorithm_official_Phase2TrackerDigi_"),
#            MinPtThreshold = cms.double(0.8)
#        )
    return process

# Execute the customization
process = apply_upc_low_pt_settings(process)
'''

# --- L1 TRACK NTUPLE MAKER START ---
process.load("L1Trigger.TrackFindingTracklet.L1TrackHitNtupleMaker_cfi")
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')

# Point the ntuple maker to your specific collections
process.L1TrackHitNtupleMaker.L1TrackInputTag = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks")
process.L1TrackHitNtupleMaker.L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted")
process.L1TrackHitNtupleMaker.MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks")
process.L1TrackHitNtupleMaker.TP_minPt = cms.double(0)
process.L1TrackHitNtupleMaker.TP_minPt = cms.double(0)
process.L1TrackHitNtupleMaker.L1Tk_minNStub = cms.int32(3)    
process.L1TrackHitNtupleMaker.TP_minNStub = cms.int32(3)
process.L1TrackHitNtupleMaker.TP_minNStubLayer = cms.int32(3)
process.L1TrackHitNtupleMaker.SaveStubs = cms.bool(True)

# Define the output file name for the ntuple
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('L1TrackHitNtuple_UPC_v4.root')
#    fileName = cms.string('/eos/cms/store/group/phys_heavyions/davidlw/L1TrackTrigger/hydjet/L1TrackHitNtuple_hydjet_v4.root')
)

# Create a path for the ntuplizer
process.ntuple_step = cms.Path(process.L1TrackHitNtupleMaker)

# Updated schedule: Only keep Tracking, Simulation, and your Ntuple
process.schedule = cms.Schedule(
#    process.L1TrackTrigger_step, # Essential: Reconstructs the tracks
#    process.StubRECO_step,       # Essential: Reconstructs the stubs
    process.L1simulation_step,   # Essential: Provides the L1 objects
    process.ntuple_step,         # Your new L1TrackHitNtupleMaker path
    process.endjob_step          # Standard cleanup
)
# --- L1 TRACK NTUPLE MAKER END ---
