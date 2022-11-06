#!/bin/bash

# ConfDB configurations to use
MASTER="/dev/CMSSW_12_4_0/HLT"      # no explicit version, take the most recent
TARGET="/dev/CMSSW_12_4_0/\$TABLE"  # no explicit version, take the most recent

TABLES="GRun HIon PIon PRef"        # $TABLE in the above variable will be expanded to these TABLES

# command-line arguments
VERBOSE=false # print extra messages to stdout
DBPROXYOPTS="" # db-proxy configuration
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v) VERBOSE=true; shift;;
    -q) VERBOSE=false; shift;;
    --dbproxy) DBPROXYOPTS="${DBPROXYOPTS} --dbproxy"; shift;;
    --dbproxyhost) DBPROXYOPTS="${DBPROXYOPTS} --dbproxyhost $2"; shift; shift;;
    --dbproxyport) DBPROXYOPTS="${DBPROXYOPTS} --dbproxyport $2"; shift; shift;;
    *) shift;;
  esac
done

# remove spurious whitespaces and tabs from DBPROXYOPTS
DBPROXYOPTS=$(echo "${DBPROXYOPTS}" | xargs)

# path to directory hosting this script
TESTDIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

# this is used for brace expansion
TABLES_=$(echo $TABLES | sed -e's/ \+/,/g')

function log() {
  $VERBOSE && echo -e "$@"
}

function getConfig() {
  local CONFIG="$1"
  local NAME="$2"
  log "  dumping HLT cffs for $NAME from $CONFIG"

  # do not use any conditions or L1 override
  hltGetConfiguration --cff --data ${CONFIG} --type ${NAME} ${DBPROXYOPTS} > HLT_${NAME}_cff.py
}

function getEventContent() {
  log "  dumping EventContent"
  local CONFIG="$1"
  local TARGET="$2"
  ${TESTDIR}/getEventContent.py ${CONFIG} ${DBPROXYOPTS} > ${TARGET}
}

function getDatasets() {
  log "  dumping Primary Dataset"
  local CONFIG="$1"
  local TARGET="$2"
  ${TESTDIR}/getDatasets.py ${CONFIG} ${DBPROXYOPTS} > ${TARGET}
}

function getConfigForOnline() {
  local CONFIG="$1"
  local NAME="$2"
# local L1T="tag[,connect]" - record is hardwired as L1GtTriggerMenuRcd

# local L1TPP="L1GtTriggerMenu_L1Menu_Collisions2012_v3_mc,sqlite_file:/afs/cern.ch/user/g/ghete/public/L1Menu/L1Menu_Collisions2012_v3/sqlFile/L1Menu_Collisions2012_v3_mc.db"
# local L1TPP="L1GtTriggerMenu_L1Menu_Collisions2012_v3_mc"
# local L1TPP="L1GtTriggerMenu_L1Menu_Collisions2015_25ns_v1_mc,sqlite_file:/afs/cern.ch/user/g/ghete/public/L1Menu/L1Menu_Collisions2015_25ns_v1/sqlFile/L1Menu_Collisions2015_25ns_v1_mc.db"
# local L1THI="L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2011_v0_mc,sqlite_file:/afs/cern.ch/user/g/ghete/public/L1Menu/L1Menu_CollisionsHeavyIons2011_v0/sqlFile/L1Menu_CollisionsHeavyIons2011_v0_mc.db"
# local L1THI="L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2011_v0_mc"
# local L1THI="L1GtTriggerMenu_L1Menu_Collisions2012_v3_mc"
# local L1THI="L1GtTriggerMenu_L1Menu_Collisions2015_25ns_v1_mc,sqlite_file:/afs/cern.ch/user/g/ghete/public/L1Menu/L1Menu_Collisions2015_25ns_v1/sqlFile/L1Menu_Collisions2015_25ns_v1_mc.db"
# local L1TPI="L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2013_v0_mc,sqlite_file:/afs/cern.ch/user/g/ghete/public/L1Menu/L1Menu_CollisionsHeavyIons2013_v0/sqlFile/L1Menu_CollisionsHeavyIons2013_v0_mc.db"
# local L1TPI="L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2013_v0_mc"
# local L1TPI="L1GtTriggerMenu_L1Menu_Collisions2012_v3_mc"
# local L1TPI="L1GtTriggerMenu_L1Menu_Collisions2015_25ns_v1_mc,sqlite_file:/afs/cern.ch/user/g/ghete/public/L1Menu/L1Menu_Collisions2015_25ns_v1/sqlFile/L1Menu_Collisions2015_25ns_v1_mc.db"

  local L1TPP1=""
  local L1TPP2=""

  log "  dumping full HLT for ${NAME} from ${CONFIG}"

  # override L1 menus
  local AUTOGT="auto:run3_hlt_${NAME}"
  if [ "${NAME}" = "Fake1" ] || [ "${NAME}" = "Fake2" ] || [ "${NAME}" = "2018" ]; then
    AUTOGT="auto:run2_hlt_${NAME}"
  elif [ "${NAME}" = "Fake" ]; then
    AUTOGT="auto:run1_hlt_${NAME}"
  fi

  hltGetConfiguration --full --data "${CONFIG}" --type "${NAME}" --unprescale --process "HLT${NAME}" --globaltag "${AUTOGT}" \
    --input "file:RelVal_Raw_${NAME}_DATA.root" ${DBPROXYOPTS} > OnLine_HLT_"${NAME}".py
}

# make sure we're using *this* working area
eval `scram runtime -sh`

# cff fragments in HLTrigger/Configuration/python/
echo "Extracting cff python dumps"
FILES=$(eval echo HLT_FULL_cff.py HLT_{$TABLES_}_cff.py HLTrigger_Datasets_{$TABLES_}_cff.py HLTrigger_EventContent_cff.py )
rm -f ${FILES}
getConfig ${MASTER} FULL
getEventContent ${MASTER} HLTrigger_EventContent_cff.py
for TABLE in ${TABLES}; do
  log "${TABLE}"
  echo "${TABLE}"
  getConfig $(eval echo ${TARGET}) ${TABLE}
  getDatasets $(eval echo ${TARGET}) HLTrigger_Datasets_${TABLE}_cff.py
done
log "Done"
log "$(ls -l ${FILES})"
mv -f ${FILES} ../python/
log

# full configs in HLTrigger/Configuration/test/
log "Extracting full configuration dumps"
echo "Extracting full configuration dumps"
FILES=$(eval echo OnLine_HLT_FULL.py OnLine_HLT_{$TABLES_}.py)
rm -f ${FILES}
getConfigForOnline ${MASTER} FULL
for TABLE in ${TABLES}; do
  log "${TABLE}"
  echo "${TABLE}"
  getConfigForOnline $(eval echo ${TARGET}) ${TABLE}
done
log "Done"
log "$(ls -l ${FILES})"
log
