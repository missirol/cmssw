#!/bin/bash

# Pass in name and status
function die {
  printf "\n%s: status %s\n" "$1" "$2"
  if [ $# -gt 2 ]; then
    printf "%s\n" "=== Log File =========="
    cat $3
    printf "%s\n" "=== End of Log File ==="
  fi
  exit $2
}

# run test job
TESTDIR="${LOCALTOP}"/src/DataFormats/Scouting/test

cmsRun "${TESTDIR}"/scoutingCollectionsIO_cfg.py -- \
  -i /store/mc/Run3Summer22DR/GluGlutoHHto2B2Tau_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/GEN-SIM-RAW/124X_mcRun3_2022_realistic_v12-v2/2550000/bbfb86f3-4073-47e3-967b-059aa6b904ad.root \
  -n 150 --skip 0 -o testDataFormatsScoutingRun3_step1.root &> testDataFormatsScoutingRun3_step1.log \
  || die "Failure running scoutingCollectionsIO_cfg.py" $? testDataFormatsScoutingRun3_step1.log

cat testDataFormatsScoutingRun3_step1.log

# compare to expected output of test job
"${TESTDIR}"/scoutingCollectionsDumper.py -v 1 -n 1 -s 137 -i testDataFormatsScoutingRun3_step1.root -k Run3Scouting \
  &> testDataFormatsScoutingRun3_step2.log \
  || die "Failure running scoutingCollectionsDumper.py" $? testDataFormatsScoutingRun3_step2.log

diff -q "${TESTDIR}"/testDataFormatsScoutingRun3_expected.log testDataFormatsScoutingRun3_step2.log \
  || die "Unexpected differences in outputs of testDataFormatsScoutingRun3 (step 2)" $?
