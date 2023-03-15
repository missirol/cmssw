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

cmsRun "${TESTDIR}"/testRun3Scouting_step1.py -- -n 100 -o testRun3Scouting_tmp.root &> log_testRun3Scouting_step1 \
  || die "Failure running testRun3Scouting_step1.py" $? log_testRun3Scouting_step1

cat log_testRun3Scouting_step1

# compare to expected output of test job
"${TESTDIR}"/testRun3Scouting_step2.py -v 1 -n 10 -i testRun3Scouting_tmp.root &> log_testRun3Scouting_step2 \
  || die "Failure running testRun3Scouting_step2.py" $? log_testRun3Scouting_step2

diff -q "${TESTDIR}"/testRun3Scouting_step2_output_expected.txt log_testRun3Scouting_step2 \
  || die "Unexpected differences in outputs of testRun3Scouting_step2.py" $?
