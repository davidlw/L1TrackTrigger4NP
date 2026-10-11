# Wrapper around the repo's v2 pset (L1Trigger/TrackFindingTracklet/test/rerunL1_trackhitntuple_cfg.py):
# only input, event count and output name change.
import os, sys
import FWCore.ParameterSet.Config as cms
base = os.path.join(os.environ["CMSSW_BASE"], "src/L1Trigger/TrackFindingTracklet/test/rerunL1_trackhitntuple_cfg.py")
exec(open(base).read())
infile = os.environ["DIAG_INPUT"]
process.source.fileNames = cms.untracked.vstring("file:" + infile)
process.maxEvents.input = int(os.environ.get("DIAG_NEVENTS", "-1"))
process.TFileService.fileName = cms.string("diag_ntuple.root")
if hasattr(process, "MessageLogger"):
    process.MessageLogger.cerr.FwkReport.reportEvery = 100

# The committed schedule leaves L1TrackTrigger_step commented out, so it re-makes the dummy stubs
# but neither re-runs the tracklet track finding nor the stub/track truth association (the
# ntuple then has no tracks and no stub truth). L1TrackTrigger = clusters+stubs, stub
# association, DTC, hybrid tracks with track association. The trigger-menu emulator
# (L1simulation_step) does not affect the tracks and is left out.
process.schedule = cms.Schedule(process.L1TrackTrigger_step, process.ntuple_step, process.endjob_step)

# Pass TrackingParticle labels into the tracklet emulation (needed by the truth trace only;
# the algorithm does not use them).
process.l1tTTTracksFromTrackletEmulation.readMoreMcTruth = True

# Prompt hybrid tracking only: the extended (displaced) emulation runs its own copy of the
# tracklet code and would write to the same truthtrace.txt.
process.L1TrackTrigger_step = cms.Path(process.TrackTriggerClustersStubs * process.TrackTriggerAssociatorClustersStubs
                                       * process.ProducerDTC * process.offlineBeamSpot
                                       * process.l1tTTTracksFromTrackletEmulation * process.TrackTriggerAssociatorTracks)
process.schedule = cms.Schedule(process.L1TrackTrigger_step, process.ntuple_step, process.endjob_step)
