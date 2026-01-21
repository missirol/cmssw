#!/bin/bash -e

if [ $# -lt 1 ]; then
  printf "\n%s\n\n" ">> ERROR: input argument missing - specify path to output directory"
  exit 1
fi

# number of events per sample
NEVT=100000

if [ $# -eq 1 ]; then
  ODIR=${1}
#  ODIR_cmsRun=$1
else
  ODIR=${1}
#  ODIR_cmsRun=${2}
fi

if [ -d ${ODIR} ]; then
  printf "%s\n" "output directory already exists: ${ODIR}"
  exit 1
fi

declare -A samplesMap

# QCD Pt-Flat
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block195ba5ab"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#195ba5ab-97d9-4005-a533-ed89974c5b1c"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block1fcb5cac"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#1fcb5cac-c028-46a8-8b9c-2e20868b1ef4"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block039d95e5"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#039d95e5-016e-4594-8364-4729521f3b5b"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block0b0179c7"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#0b0179c7-c911-48b2-bf17-d993c3cfef91"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block0cb27300"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#0cb27300-017c-4c9a-9526-e8979eb8b3d2"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block0cebc851"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#0cebc851-7f46-40b8-a070-06443543caf4"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block0752d210"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#0752d210-0412-4a5c-a119-b5cc9531eeab"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block0a22fa4f"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#0a22fa4f-e2fe-45ba-9c4d-c3a7d47f8217"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block1615a814"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#1615a814-831c-4536-8c80-e811ce70ef57"
samplesMap["Run3Winter25_QCD_PtFlat15to7000_13p6TeV_FlatPU0to120_block0d60974e"]="/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/Run3Winter25Digi-FlatPU0to120_142X_mcRun3_2025_realistic_v9-v4/GEN-SIM-RAW#0d60974e-d067-4b29-8bd0-5453b264ea01"

# options (JobFlavour and AccountingGroup)
opts=""
if [[ ${HOSTNAME} == lxplus* ]]; then
  opts+="--JobFlavour espresso"
fi

COMMON_OPTS=" --filein tmp.root"
COMMON_OPTS+=" --mc --conditions auto:phase1_2025_realistic --geometry DB:Extended"
COMMON_OPTS+=" --scenario pp --era Run3_2025"
COMMON_OPTS+=" --datatier NANOAOD --eventcontent NANOAOD"
COMMON_OPTS+=" --nThreads 8 --nStreams 0"
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

rm -rf "${JOB_LABEL}"_cfg.py __pycache__

for sampleKey in ${!samplesMap[@]}; do
  sampleName=${samplesMap[${sampleKey}]}

  # number of events per sample
  numEvents=${NEVT}

  bdriver -c "${JOB_LABEL}"_cfg_dump.py --customize-cfg -m ${numEvents} -n 1000 --cpus 8 --mem 1000 --time 600 ${opts} \
    -d ${sampleName} -p 0 -o ${ODIR}/${sampleKey}
done
unset sampleKey numEvents sampleName

unset opts samplesMap NEVT ODIR #ODIR_cmsRun
