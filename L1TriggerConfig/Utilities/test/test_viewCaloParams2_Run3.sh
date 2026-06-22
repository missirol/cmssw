#!/bin/bash

if [ -z "${SCRAM_TEST_PATH}" ]; then
  SCRAM_TEST_PATH=$(dirname -- "${BASH_SOURCE[0]}")
  set -e
fi

exec_cfi () {
  printf "%s\n" "Config: ${1}"
  cmsRun "${SCRAM_TEST_PATH}"/viewCaloParams2.py \
    -c L1Trigger.L1TCalorimeter."${1}" &> l1tCaloParams_"${1}".txt
}

exec_mc () {
  printf "%s\n" "MC GlobalTag: ${1}"
  cmsRun "${SCRAM_TEST_PATH}"/viewCaloParams2.py \
    -g "${1}" &> l1tCaloParams_"${1}".txt
}

exec_data () {
  tmpGlobalTag=160X_dataRun3_HLT_v1
  printf "%s\n" "Data GlobalTag: ${tmpGlobalTag} (run-${1}, label: \"${2}\")"
  cmsRun "${SCRAM_TEST_PATH}"/viewCaloParams2.py \
    -g "${tmpGlobalTag}" -r "${1}" &> l1tCaloParams_run"${1}"_"${2}".txt
}

# 2023 pp
exec_cfi caloParams_2023_v0_2_cfi
exec_mc 130X_mcRun3_2023_realistic_v15
exec_data 370293 2023_pp

# 2023 PbPb
exec_cfi caloParamsHI_2023_v0_4_3_cfi
exec_mc 132X_mcRun3_2023_realistic_HI_v10
exec_data 375790 2023_PbPb

# 2024 pp
exec_cfi  caloParams_2024_v0_3_cfi
exec_mc   140X_mcRun3_2024_realistic_v26
exec_data 386025 2024_pp-v1
exec_data 386924 2024_pp-v2

# 2024 PbPb
exec_cfi  caloParamsHI_2024_v0_2_cfi
exec_mc   141X_mcRun3_2024_realistic_HI_v17
exec_data 388750 2024_PbPb

# 2025 pO
exec_mc   150X_mcRun3_2025_forpO_realistic_v9
exec_data 394004 2025_pO

# 2025 OO
exec_cfi  caloParamsOO_2025_v0_0_cfi
exec_mc   150X_mcRun3_2025_forOO_realistic_v9
exec_data 394217 2025_OO

# 2025 NeNe
exec_mc   150X_mcRun3_2025_forNeNe_realistic_v9
exec_data 394272 2025_NeNe

# 2025 pp
exec_cfi  caloParams_2025_v0_3_cfi
exec_mc   150X_mcRun3_2025_realistic_v14
exec_data 394959 2025_pp-v1
exec_data 398860 2025_pp-v2

# 2025 PbPb
exec_cfi  caloParamsHI_2025_v0_0_cfi
exec_mc   151X_mcRun3_2025_realistic_HI_v5
exec_data 400391 2025_PbPb

# 2026 pp
exec_mc 160X_mcRun3_2026_realistic_v6
exec_data 403937 2026_pp

# 2026 PbPb
exec_data 404925 2026_PbPb
