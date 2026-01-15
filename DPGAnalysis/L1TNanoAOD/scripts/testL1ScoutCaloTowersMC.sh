#!/bin/bash

INPUT_FILE=/store/mc/Run3Winter25Digi/SinglePion_E-50_Eta-0to3-pythia8-gun/GEN-SIM-RAW/NoPU_142X_mcRun3_2025BOY_realistic_Candidate_2024_11_13_17_21_33-v2/2530000/0515ebc7-845c-4d78-8c34-f5ce68d68164.root

COMMON_OPTS=" --filein ${INPUT_FILE}"
COMMON_OPTS+=" --mc --conditions auto:phase1_2025_realistic --geometry DB:Extended"
COMMON_OPTS+=" --scenario pp --era Run3_2025"
COMMON_OPTS+=" --datatier NANOAOD --eventcontent NANOAOD"
COMMON_OPTS+=" --nThreads 1 --nStreams 0"
COMMON_OPTS+=" --no_exec"

JOB_LABEL=nanoL1TCustom
cmsDriver.py "${JOB_LABEL}" --process "${JOB_LABEL^^}" ${COMMON_OPTS} \
  --python_filename "${JOB_LABEL}"_cfg.py --fileout file:"${JOB_LABEL}"_out.root \
  -s RAW2DIGI,NANO:@GENLite+@L1ScoutCaloTowersMC \
  -n 10

cat <<@EOF >> "${JOB_LABEL}"_cfg.py
process.NANOAODoutput.saveTriggerResults = cms.untracked.bool(False)
@EOF

edmConfigDump --prune "${JOB_LABEL}"_cfg.py > "${JOB_LABEL}"_cfg_dump.py

cmsRun "${JOB_LABEL}"_cfg_dump.py 2>&1 | tee "${JOB_LABEL}"_cfg_dump.log

rm -rf __pycache__
