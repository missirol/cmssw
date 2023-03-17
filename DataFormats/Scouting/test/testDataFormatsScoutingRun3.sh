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

cmsRun "${TESTDIR}"/testDataFormatsScoutingRun3_step1.py -- \
  -i /store/mc/Run3Summer22DR/GluGlutoHHto2B2Tau_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/GEN-SIM-RAW/124X_mcRun3_2022_realistic_v12-v2/2550000/bbfb86f3-4073-47e3-967b-059aa6b904ad.root \
  -n 150 -o testDataFormatsScoutingRun3_tmp.root &> log_testDataFormatsScoutingRun3_step1 \
  || die "Failure running testDataFormatsScoutingRun3_step1.py" $? log_testDataFormatsScoutingRun3_step1

cat log_testDataFormatsScoutingRun3_step1

# compare to expected output of test job
"${TESTDIR}"/testDataFormatsScoutingRun3_step2.py -v 1 -n 1 -s 137 -i testDataFormatsScoutingRun3_tmp.root \
  &> log_testDataFormatsScoutingRun3_step2 \
  || die "Failure running testDataFormatsScoutingRun3_step2.py" $? log_testDataFormatsScoutingRun3_step2

diff -q "${TESTDIR}"/testDataFormatsScoutingRun3_step2_output_expected.txt log_testDataFormatsScoutingRun3_step2 \
  || die "Unexpected differences in outputs of testDataFormatsScoutingRun3_step2.py" $?
