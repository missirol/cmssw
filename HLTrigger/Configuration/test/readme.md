The `test/` directory of the `HLTrigger/Configuration` package contains

 - scripts to copy HLT menus from the `ConfDB` database into CMSSW,
   as both `cff` fragments (loadable via `cmsDriver.py`) and standalone `cfg` configurations (usable with `cmsRun`);

 - scripts to run tests with these HLT menus
   (a version of these tests runs in CMSSW integration buils as the so-called "HLT-Validation" tests);

 - a unit test to verify the availability of the EDM input files required by the main set of HLT tests running in CMSSW
   (see `test_AccessToEDMInputsOfHLTTests` below).

---

Unit test: `test_AccessToEDMInputsOfHLTTests`
---

This unit test executes `cmsRun` jobs to verify the availability
of the EDM files listed in `HLTrigger/Configuration/test/testAccessToEDMInputsOfHLTTests_filelist.txt`.

To run the unit test via `scram`, execute
```sh
scram b runtests_test_AccessToEDMInputsOfHLTTests
```
To run the unit test locally, execute
```sh
LOCALTOP="${CMSSW_BASE}" "${CMSSW_BASE}"/src/HLTrigger/Configuration/test/testAccessToEDMInputsOfHLTTests.sh
```

The unit test does not modify the content of the file `testAccessToEDMInputsOfHLTTests_filelist.txt`.
The latter can be updated by manually executing the script `testAccessToEDMInputsOfHLTTests_update_filelist.sh`.
The file `testAccessToEDMInputsOfHLTTests_filelist.txt` lists
the Logical File Name (LFN) of the EDM files used in HLT tests
for the main CMSSW development branches (name format: `CMSSW_[0-9]*_[0-9]*_X`).
The list includes only EDM files which are available either in the IB-EOS cache at the CERN T2, or at any T2/T3 site.
