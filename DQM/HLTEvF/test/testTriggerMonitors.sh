#!/bin/bash

# Pass in name and status
function die {
  echo $1: status $2
  echo === Log file ===
  cat ${3:-/dev/null}
  echo === End log file ===
  exit $2
}

# run test job
TESTDIR="${LOCALTOP}"/src/DQM/HLTEvF/test

cmsRun "${TESTDIR}"/testTriggerMonitors_dqm_cfg.py &> log_testTriggerMonitors_dqm \
  || die "Failure running testTriggerMonitors_dqm_cfg.py" $? log_testTriggerMonitors_dqm

cmsRun "${TESTDIR}"/testTriggerMonitors_harvesting_cfg.py &> log_testTriggerMonitors_harvesting \
  || die "Failure running testTriggerMonitors_harvesting_cfg.py" $? log_testTriggerMonitors_harvesting
