# ---------------------------------------------------------------------------
# offlineTrackNtuple_cfg.py
#
# Offline tracks matched to TrackingParticles -> efficiency, fake rate, duplicates.
#
#     cmsRun offlineTrackNtuple_cfg.py
#
# INPUT: step3.root only. It must have been written with the extra keeps in
# simulations/commands_cmsDriver -- the association is hit-based and needs
# generalTracks + its trackExtras + rec hits, plus siPixelClusters and
# siPhase2Clusters. The TrackingParticles and DigiSimLinks come along from step 2
# via the blanket 'keep *_*_*_HLT'. Check before a long run:
#
#     edmDumpEventContent step3.root | grep -E "generalTracks|Cluster|TrackingParticle"
#
# If the keeps are missing the job still runs and every trk_* branch is EMPTY.
#
# This job does NOT re-emulate L1 and contains no L1 code; the L1 ntuple is made
# separately by rerunL1_trackhitntuple_cfg.py. Join the two by event if needed.
# ---------------------------------------------------------------------------

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('OFFLINETRK', Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun4_realistic_v4', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step3.root'),
)

# --- offline track <-> TrackingParticle association ------------------------
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

# step3.root carries generalTracks, siPixelClusters and siPhase2Clusters TWICE:
# from step 2's HLT step (process HLT -- the trigger's own online tracking) and
# from step 3 (process RECO -- the offline reconstruction). Pin to RECO. An
# unqualified tag resolves to RECO today, but relying on that would let a silent
# switch to HLT tracks pass unnoticed: the numbers would still look sane.
OFFLINE = "RECO"
process.tpClusterProducer.pixelClusterSrc     = cms.InputTag("siPixelClusters",  "", OFFLINE)
process.tpClusterProducer.phase2OTClusterSrc  = cms.InputTag("siPhase2Clusters", "", OFFLINE)
process.tpClusterProducer.pixelSimLinkSrc     = cms.InputTag("simSiPixelDigis", "Pixel",   "HLT")
process.tpClusterProducer.phase2OTSimLinkSrc  = cms.InputTag("simSiPixelDigis", "Tracker", "HLT")
process.tpClusterProducer.trackingParticleSrc = cms.InputTag("mix", "MergedTrackTruth", "HLT")

process.trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("generalTracks", "", OFFLINE)
process.trackingParticleRecoTrackAsssociation.label_tp = cms.InputTag("mix", "MergedTrackTruth", "HLT")

# --- the ntuplizer ---------------------------------------------------------
process.OfflineTrackNtupleMaker = cms.EDAnalyzer('OfflineTrackNtupleMaker',
    TrackInputTag            = cms.InputTag("generalTracks", "", OFFLINE),
    AssociatorInputTag       = cms.InputTag("trackingParticleRecoTrackAsssociation"),
    BeamSpotInputTag         = cms.InputTag("offlineBeamSpot", "", OFFLINE),
    TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth", "HLT"),

    # Loose on purpose. The denominator is set upstream by ptMinTP in the mixing
    # module at step 2, not here, and carries NO stub/hit requirement -- so the
    # efficiency denominator stays independent of the reconstruction.
    TP_minPt       = cms.double(0.0),
    TP_maxEta      = cms.double(2.5),
    TP_maxZ0       = cms.double(30.0),
    TP_minNHits    = cms.int32(0),
    TP_onlyCharged = cms.bool(True),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('OfflineTrackNtuple.root'))

process.p = cms.Path(
    process.tpClusterProducer
    * process.quickTrackAssociatorByHits
    * process.trackingParticleRecoTrackAsssociation
    * process.OfflineTrackNtupleMaker
)
process.schedule = cms.Schedule(process.p)
