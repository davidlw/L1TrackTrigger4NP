## L1TrackTrigger4NP Code repository for developing low-pT L1 track trigger for nuclear collisions at Phase-2 CMS detector for HL-LHC

cd /yourepath/ # usually a CMSSW src directory

git clone https://github.com/davidlw/L1TrackTrigger4NP

## Getting started without lxplus (analysis only)

The macros under `analyses/` are plain ROOT and do not need CMSSW or lxplus.
Only ntuple *production* (the `cmsRun` steps below) does. To compute efficiency,
fake rate, duplicate rate or the trigger efficiency on a laptop:

1. Install ROOT (`brew install root`, or `conda install -c conda-forge root`).
2. `git clone https://github.com/davidlw/L1TrackTrigger4NP`
3. Get a few ntuple files from someone with EOS access -- one L1 file
   (`L1TrackHitNtuple_*.root`) is ~12 MB, so ten of them per sample is plenty
   for a first look. Put them in any directory; the file names do not matter.
4. Run a fill macro on that directory (or on one file), then the matching
   `plot_*` macro on its output:

```bash
cd analyses/offline_tracking_performance/scripts
root -l -b -q 'offline_perf.C("/my/ntuples/QED_mumu","../output/offperf_qed_mumu.root",0,false,true,4,1)'
root -l -b -q 'plot_offline_perf.C("../output/offperf_qed_mumu.root")'
```

```bash
cd analyses/L1_tracking_performance/output
root -l -b -q '../scripts/compare_eff.C(0,true,"eff_qed_mumu.root","/my/ntuples/DefaultStub/QED_mumu","/my/ntuples/DummyStub/QED_mumu")'
root -l -b -q '../scripts/plot_eff.C("eff_qed_mumu.root")'
```

Each `analyses/*/README.md` lists every macro and argument. The fill macros
print which files they found, so a wrong path shows up immediately.

## Simulations
### Available generator fragments are kept here: https://github.com/davidlw/genproductions. 
### For UPC2024 as an example: https://github.com/davidlw/genproductions/tree/UPC2024/genfragments/PbPb_5p36TeV/Starlight

cmsrel CMSSW_14_0_6

cd CMSSW_14_0_6/src 

cmsenv

mkdir -p Configuration/GenProduction/python # put fragments you plan to use here

### look at L1TrackTrigger4NP/simulations/commands_cmsDriver for cmsDriver commands to generate configurations for various steps
### step 1 for GEN-SIM; step 2 for DIGI-RAW - this step produces the output needed for L1 track trigger emulation; step 3: offline reco
### example crab config files can also be found in L1TrackTrigger4NP/simulations

## IMPORTANT for low-pT studies: the truth pT cut lives in step 2

The TrackingParticle collection is filtered in the mixing module at DIGI time, not
in the ntuplizer. The customisation `SimGeneral/MixingModule/customiseStoredTPConfig.higherPtTP`
sets `process.mix.digitizers.mergedtruth.select.ptMinTP = 1.0` GeV (CMSSW default: 0.1).

Samples produced with it contain almost no truth particles below 1 GeV, so tracking
efficiency below 1 GeV cannot be measured from them at all - the denominator is
missing. This is upstream of both the stub definition and the ntuplizer, so neither
`TP_minPt` nor the dummy-stub configuration can recover it.

`commands_cmsDriver` now sets `ptMinTP` explicitly and keeps the old command
commented out for reproducing the existing samples. Regenerating from step 2 is
required to study the sub-GeV region.

## Producing L1 track ntuple with clusters information:

cmsrel CMSSW_15_1_0_patch3

cd CMSSW_15_1_0_patch3/src

cmsenv

git cms-addpkg L1Trigger/TrackFindingTracklet

git cms-addpkg L1Trigger/TrackTrigger/

git cms-addpkg L1Trigger/TrackerDTC/

git cms-addpkg L1Trigger/TrackFindingTMTT/

git cms-addpkg SimTracker/TrackTriggerAssociation/

scram b -j4

cp L1TrackTrigger4NP/L1Trigger/TrackFindingTracklet/python/* L1Trigger/TrackFindingTracklet/python/

cp L1TrackTrigger4NP/L1Trigger/TrackFindingTracklet/test/* L1Trigger/TrackFindingTracklet/test/ 

cp L1TrackTrigger4NP/L1Trigger/TrackFindingTracklet/interface/* L1Trigger/TrackFindingTracklet/interface/

cp L1TrackTrigger4NP/L1Trigger/TrackFindingTracklet/plugins/* L1Trigger/TrackFindingTracklet/plugins/

cp L1TrackTrigger4NP/L1Trigger/TrackFindingTracklet/src/* L1Trigger/TrackFindingTracklet/src/

cp L1TrackTrigger4NP/L1Trigger/TrackTrigger/plugins/* L1Trigger/TrackTrigger/plugins

cp L1TrackTrigger4NP/L1Trigger/TrackTrigger/python/* L1Trigger/TrackTrigger/python/

cp L1TrackTrigger4NP/L1Trigger/TrackFindingTMTT/src/* L1Trigger/TrackFindingTMTT/src/

cp L1TrackTrigger4NP/L1Trigger/TrackerDTC/src/* L1Trigger/TrackerDTC/src/

cp L1TrackTrigger4NP/SimTracker/TrackTriggerAssociation/plugins/* SimTracker/TrackTriggerAssociation/plugins/

cmsenv

scram b -j4

cd L1TrackTrigger4NP/configurations

cmsRun rerunL1_trackhitntuple_cfg.py 
