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
  hlt8677a2aa_baseline
  hlt8677a2aa_CaloRecHits
  hlt8677a2aa_XcalRecHits
)

for hltLabel in "${hltLabels[@]}"; do

  merge_files "${hltLabel}"

done
unset hltLabel
