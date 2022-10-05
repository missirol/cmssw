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
TESTDIR="${LOCALTOP}"/src/L1Trigger/L1TGlobal/test

cmsRun "${TESTDIR}"/testL1TGlobalProducer_cfg.py &> log_testL1TGlobalProducer \
 || die "Failure running testL1TGlobalProducer_cfg.py" $? log_testL1TGlobalProducer

# expected PathSummary of test job
cat <<@EOF > log_testL1TGlobalProducer_expected
==================  L1 Trigger Report  =====================================================================

 L1T menu Name   : L1Menu_Collisions2022_FracPrescale_Test
 L1T menu Version: 0.10
 L1T menu Comment: Test menu with Mu5, EG, jet, ZB, and MB seeds. 

    Bit                  Algorithm Name                  Init    PScd  Final   PS Factor     Num Bx Masked
============================================================================================================
       0                               L1_SingleMu3       489    323    323       1.5          0
       1                               L1_SingleEG3       984    654    654       1.5          0
       2                              L1_SingleJet8      1000    660    660       1.5          0
       3                                L1_ZeroBias      1000    660    660       1.5          0
       4                          L1_MinimumBiasHF0      1000    660    660       1.5          0
                                                      Final OR Count = 893
@EOF

# compare to expected output of test job
sed -n '/L1 Trigger Report  =/,/Final OR Count =/p' log_testL1TGlobalProducer \
 | diff log_testL1TGlobalProducer_expected - \
 || die "differences in expected log report" $?
