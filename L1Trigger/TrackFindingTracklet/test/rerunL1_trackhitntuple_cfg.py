
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
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/phys_heavyions/davidlw/StarLightJpsiPhase2_PrivateMC/Step2_DIGI_RAW_CMSSW_14_0_6/260325_123236/0000/step2_1.root'
    ),
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

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Path and EndPath definitions
process.StubRECO_step = cms.Path(process.TrackTriggerClustersStubs)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger) # RESTORED: Reconstructs L1 tracks from stubs
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Automatic addition of the customisation functions
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000 
process = customise_aging_1000(process)

from Configuration.DataProcessing.Utils import addMonitoring 
process = addMonitoring(process)

from L1Trigger.Configuration.customisePhase2TTOn110 import customisePhase2TTOn110 
process = customisePhase2TTOn110(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)


# --- L1 TRACK NTUPLE MAKER START ---
process.load("L1Trigger.TrackFindingTracklet.L1TrackHitNtupleMaker_cfi")
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')

# Point the ntuple maker to your specific collections
process.L1TrackHitNtupleMaker.L1TrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks")
process.L1TrackHitNtupleMaker.L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted")
process.L1TrackHitNtupleMaker.MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks")

# Cluster Input tags for your custom GNN extraction setup
process.L1TrackHitNtupleMaker.L1ClusterInputTag = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive")
process.L1TrackHitNtupleMaker.phase2OTClusters = cms.InputTag("siPhase2Clusters")

process.L1TrackHitNtupleMaker.TP_minPt = cms.double(0)
process.L1TrackHitNtupleMaker.L1Tk_minNStub = cms.int32(3)    
process.L1TrackHitNtupleMaker.TP_minNStub = cms.int32(3)
process.L1TrackHitNtupleMaker.TP_minNStubLayer = cms.int32(3)
process.L1TrackHitNtupleMaker.SaveStubs = cms.bool(True)

# Define the output file name for the ntuple
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('L1TrackHitNtuple_UPC_v4.root')
)

# Create a path for the ntuplizer
process.ntuple_step = cms.Path(process.L1TrackHitNtupleMaker)

# --- ADD OUTPUT MODULE FOR CLUSTERS & STUBS START ---
process.RAWOUTF = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('TrackTrigger_ClustersStubs_Output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *', 
        'keep *_TTClustersFromPhase2TrackerDigis_*_*', 
        'keep *_TTStubsFromPhase2TrackerDigis_*_*',    
        'keep *_mix_MergedTrackTruth_*'                
    )
)

process.output_step = cms.EndPath(process.RAWOUTF)
# --- ADD OUTPUT MODULE FOR CLUSTERS & STUBS END ---

# Updated schedule: Running tracking clusters/stubs, tracks, ntuplizer, and output steps
process.schedule = cms.Schedule(
    process.StubRECO_step,       # Step 1: Clusters & Stubs creation
#    process.L1TrackTrigger_step, # Step 2: RESTORED Tracklet/KF tracking logic loop
    process.L1simulation_step,   # Step 3: Trigger primitive updates
    process.ntuple_step,         # Step 4: Your analyzer tree processing
#    process.output_step,         # Step 5: Standalone EDM formatting dump
    process.endjob_step          # Step 6: Standard cleanup
)
# --- L1 TRACK NTUPLE MAKER END ---
