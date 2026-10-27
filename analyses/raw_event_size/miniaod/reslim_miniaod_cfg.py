"""Re-write existing RECO+PAT (step-3) files with exactly the standard Phase-2
miniAOD event content and output settings, to measure the miniAOD event size.

The step-3 files were written with extra keeps (all HLT/SIM products, digis,
tracks, clusters), so their own size is not the miniAOD size. Nothing is
re-reconstructed: products are only selected and re-compressed, with the same
event content (EventContent_cff) and output customisation
(miniAOD_tools.miniAOD_customizeOutput) cmsDriver uses for the MINIAOD tiers.

Two outputs:
  <out>_sim.root  : MINIAODSIM (includes MC truth: gen particles, gen jets, ...)
  <out>_data.root : MINIAOD    (data-like content, the relevant size for real data)

Usage: cmsRun reslim_miniaod_cfg.py inputFiles=file:a.root[,...] outputFile=out.root maxEvents=N
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
from Configuration.EventContent.EventContent_cff import MINIAODSIMEventContent, MINIAODEventContent
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeOutput

opts = VarParsing("analysis")
opts.parseArguments()
base = opts.outputFile.replace(".root", "")

process = cms.Process("SLIM", Phase2C17I13M9)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(opts.inputFiles))

for label, content in (("sim", MINIAODSIMEventContent), ("data", MINIAODEventContent)):
    out = cms.OutputModule("PoolOutputModule", content.clone(),
                           fileName=cms.untracked.string(f"{base}_{label}.root"))
    miniAOD_customizeOutput(out)
    setattr(process, f"out_{label}", out)

process.e = cms.EndPath(process.out_sim + process.out_data)
