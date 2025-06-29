#!/bin/bash -e

merge_files() {
  hltLabel="${1}"
  outFile="${hltLabel}".root
  [ ! -f "${outFile}" ] || return
  ls "${hltLabel}"_run*.root > "${hltLabel}".txt
  edmCopyPickMerge filePrepend=file: outputFile="${outFile}" inputFiles_load="${hltLabel}".txt
  rm -f "${hltLabel}".txt
}

hltLabels=(
  hltCaloTowers250629_baseline
  hltCaloTowers250629_CaloTowers
  hltCaloTowers250629_CaloTowers_Egt1p0
  hltCaloTowers250629_CaloTowers_Egt1p0_AbsEtaLt3p0
  hltCaloTowers250629_PFRecHits
  hltCaloTowers250629_PFRecHits_Egt1p0
)

for hltLabel in "${hltLabels[@]}"; do

  merge_files "${hltLabel}"

done
unset hltLabel
