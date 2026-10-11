# usage: run_one.sh <step2 file index> <nevents>
i=$1; n=$2
D=${DIAG_INDIR:-/eos/cms/store/group/phys_heavyions/davidlw/StarLightQEDMuMuPhase2_PrivateMC/Step2_DIGI_RAW_CMSSW_14_0_6_v4_tpminpt0p3/260815_083524/0000}
J=/tmp/davidlw_l1diag/runs/${DIAG_TAG:-job}$i; mkdir -p $J && cd $J
source /cvmfs/cms.cern.ch/cmsset_default.sh; export SCRAM_ARCH=el8_amd64_gcc12
cd /tmp/davidlw_l1diag/CMSSW_15_1_0_patch3/src && eval $(scram runtime -sh) && cd - > /dev/null
export DIAG_INPUT=$D/step2_$i.root DIAG_NEVENTS=$n
cmsRun ../diag_cfg.py > cmsrun.log 2>&1; echo "exit $?" >> cmsrun.log
