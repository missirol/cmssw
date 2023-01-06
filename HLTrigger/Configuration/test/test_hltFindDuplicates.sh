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

if [ -z "${SCRAM_TEST_PATH}" ]; then
  printf "\n%s\n\n" "ERROR -- environment variable SCRAM_TEST_PATH not defined"
  exit 1
fi

# run test job
TEST_MENUS=(
#  Fake
#  Fake1
#  Fake2
#  FULL
#  GRun
#  HIon
  PRef
#  PIon
)

for TEST_MENU in "${TEST_MENUS[@]}"; do
  hltFindDuplicates "${SCRAM_TEST_PATH}"/OnLine_HLT_"${TEST_MENU}".py -x realData=0 globalTag=@ \
    -o test_hltFindDuplicates_"${TEST_MENU}"_output &> test_hltFindDuplicates_"${TEST_MENU}"_log \
    || die "Failure running hltFindDuplicates (menu: ${TEST_MENU})" $? test_hltFindDuplicates_"${TEST_MENU}"_log
done
