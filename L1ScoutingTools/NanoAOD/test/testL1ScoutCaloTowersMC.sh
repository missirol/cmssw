#!/bin/bash

INPUT_FILE=/store/mc/Run3Winter25Digi/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/2520006/18ee11ae-e5e7-4091-a73b-bd33e889a2d5.root

cd $(dirname -- "${BASH_SOURCE[0]}")

JOB_LABEL=testL1ScoutCaloTowersMC

COMMON_OPTS=" --filein ${INPUT_FILE}"
COMMON_OPTS+=" --mc --conditions auto:phase1_2025_realistic --geometry DB:Extended"
COMMON_OPTS+=" --scenario pp --era Run3_2025"
COMMON_OPTS+=" --datatier NANOAOD --eventcontent NANOAOD"
COMMON_OPTS+=" --nThreads 1 --nStreams 0"
COMMON_OPTS+=" --no_exec"

cmsDriver.py "${JOB_LABEL}" --process "${JOB_LABEL^^}" ${COMMON_OPTS} \
  --python_filename "${JOB_LABEL}"_cfg.py --fileout file:"${JOB_LABEL}"_out.root \
  -s RAW2DIGI,NANO:@GEN+@L1ScoutCaloTowersMC \
  -n 10

cat <<@EOF >> "${JOB_LABEL}"_cfg.py
process.NANOAODoutput.saveTriggerResults = cms.untracked.bool(False)
@EOF

# --prune
edmConfigDump \
 "${JOB_LABEL}"_cfg.py > "${JOB_LABEL}"_cfg_dump.py

rm -rf "${JOB_LABEL}"_cfg.py

#cmsRun "${JOB_LABEL}"_cfg_dump.py 2>&1 | tee "${JOB_LABEL}"_cfg_dump.log

rm -rf __pycache__
