# L1TrackTrigger4NP Code repository for developing low-pT L1 track trigger for nuclear collisions at Phase-2 CMS detector for HL-LHC
cd /yourepath/ # usually a CMSSW src directory
git clone https://github.com/davidlw/L1TrackTrigger4NP

# simulations
# Available generator fragments are kept here: https://github.com/davidlw/genproductions. 
# For UPC2024 as an example: https://github.com/davidlw/genproductions/tree/UPC2024/genfragments/PbPb_5p36TeV/Starlight
cmsrel CMSSW_14_0_6
cd CMSSW_14_0_6/src 
cmsenv
mkdir -p Configuration/GenProduction/python # put fragments you plan to use here
# look at L1TrackTrigger4NP/simulations/commands_cmsDriver for cmsDriver commands to generate configurations for various steps
# step 1 for GEN-SIM; step 2 for DIGI-RAW - this step produces the output needed for L1 track trigger emulation; step 3: offline reco
# example crab config files can also be found in L1TrackTrigger4NP/simulations

#Producing L1 track ntuple with clusters information:
cmsrel CMSSW_15_1_0_patch3
cd CMSSW_15_1_0_patch3/src
cmsenv
git cms-addpkg L1Trigger/TrackFindingTracklet
cp /yourpath/L1TrackTrigger4NP/L1Trigger/TrackFindingTracklet/python/*.py L1Trigger/TrackFindingTracklet/python/
cp /yourpath/L1TrackTrigger4NP/L1Trigger/TrackFindingTracklet/test/* L1Trigger/TrackFindingTracklet/test/ 
cmsenv
scram b -j4
cd /yourpath/L1TrackTrigger4NP/configurations
cmsRun rerunL1_trackhitntuple_cfg.py # require DIGI-RAW files as inputs 
