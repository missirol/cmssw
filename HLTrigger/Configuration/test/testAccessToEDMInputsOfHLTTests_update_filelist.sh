#!/bin/bash -e

# This script updates the file "${CMSSW_BASE}"/src/HLTrigger/Configuration/test/testAccessToEDMInputsOfHLTTests_filelist.txt
# with the list of EDM files potentially used by HLT tests in the main release cycles of CMSSW (i.e. branches named CMSSW_\d_\d_X).

# path to output file
outputFile="${CMSSW_BASE}"/src/HLTrigger/Configuration/test/testAccessToEDMInputsOfHLTTests_filelist.txt

# path to directory hosting this script
TESTDIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

# ensure that directory hosting this script corresponds to ${CMSSW_BASE}/src/HLTrigger/Configuration/test
if [ "${TESTDIR}" != "${CMSSW_BASE}"/src/HLTrigger/Configuration/test ]; then
  printf "\n%s\n" "ERROR -- the directory hosting testAccessToHLTTestInputs.sh [1] does not correspond to \${CMSSW_BASE}/src/HLTrigger/Configuration/test [2]"
  printf "%s\n"   "         [1] ${TESTDIR}"
  printf "%s\n\n" "         [2] ${CMSSW_BASE}/src/HLTrigger/Configuration/test"
  exit 1
fi

# files in CMSSW using EDM inputs for HLT tests
cmsswFiles=(
  HLTrigger/Configuration/test/cmsDriver.csh
  Configuration/HLT/python/addOnTestsHLT.py
  Utilities/ReleaseScripts/scripts/addOnTests.py
)

# list of CMSSW branches to be checked
# official-cmssw is the default name of the remote corresponding to the central CMSSW repository
cmsswBranches=$(git branch -a | grep 'remotes/official-cmssw/CMSSW_[0-9]*_[0-9]*_X$')
cmsswBranches+=("HEAD") # add HEAD to include updates committed locally

# create 1st temporary file (list of EDM input files incl. duplicates)
TMPFILE1=$(mktemp)

# grep from base directory
cd "${CMSSW_BASE}"/src

# loop over CMSSW branches to be grep-d
for cmsswBranch in "${cmsswBranches[@]}"; do
  git grep -h "[='\" ]/store/.*.root" ${cmsswBranch} -- ${cmsswFiles[*]} |
    sed 's|=/store/| /store/|g' | sed "s|'| |g" | sed 's|"| |g' | \
    awk '{ for(i=1;i<=NF;i++) if ($i ~ /\/store\/.*.root/) print $i }' >> "${TMPFILE1}"
done; unset cmsswBranch

# create 2nd temporary file with list of available files (without duplicates)
TMPFILE2=$(mktemp)

# a file is considered as available if present in the ibeos cache (CERN T2), or at any T2/T3 site
for inputFile in $(cat "${TMPFILE1}" | sort -u); do
  if [ $(ls /eos/cms/store/user/cmsbuild/"${inputFile}" 2> /dev/null | wc -l) -eq 0 ]; then
    if [ $(dasgoclient -query "site file=${inputFile}" | grep -E '(T2|T3)_' | wc -l) -eq 0 ]; then
      printf "%s\n" "File not available: ${inputFile}"
      continue
    fi
  fi
  echo "${inputFile}" >> "${TMPFILE2}"
done
unset inputFile

# create/update output file
cat "${TMPFILE2}" > "${outputFile}"

# return to test/ directory
cd "${TESTDIR}"
