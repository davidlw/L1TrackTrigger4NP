import FWCore.ParameterSet.Config as cms

# Configuration for Ntuple maker for analyzing L1 track performance.
# The parameters specified here are suitable for Hybrid prompt track collections

L1TrackHitNtupleMaker = cms.EDAnalyzer('L1TrackHitNtupleMaker',
       # MyProcess is the (unsigned) PDGID corresponding to the process which is run
       # e.g. single electron/positron = 11
       #      single pion+/pion- = 211
       #      single muon+/muon- = 13
       #      pions in jets = 6
       #      taus = 15
       #      all TPs = 1 (pp collisions)
       MyProcess = cms.int32(1),
       DebugMode = cms.bool(False),      # printout lots of debug statements
       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
       SaveStubs = cms.bool(True),       # CHANGED TO TRUE: Crucial to keep open for saving your custom GNN cluster vectors!
       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 tracking?
       L1Tk_minNStub = cms.int32(3),     # CHANGED TO 3: To accept soft tracks for UPC
       TP_minNStub = cms.int32(3),       # CHANGED TO 3: Denominator tracking matching boundary
       TP_minNStubLayer = cms.int32(3),  # CHANGED TO 3: 
       TP_minPt = cms.double(0.0),       # CHANGED TO 0.0: Essential to open the gate for low-pt tracks!
       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
       
       L1TrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"), # TTTrack input
       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), # MCTruth input
       
       # === UPDATED AND NEW DATA COLLECTIONS ===
       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"), 
       L1ClusterInputTag = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive"), # ADDED: Default to full raw data inclusive pool for GNN
       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
       # ========================================

       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
       TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
       
       # tracking in jets (--> requires AK4 genjet collection present!)
       TrackingInJets = cms.bool(False),
       GenJetInputTag = cms.InputTag("ak4GenJets", ""),

       phase2OTDigis = cms.InputTag("mix", "Tracker"),
       phase2OTClusters = cms.InputTag("siPhase2Clusters")
)
