#!/bin/bash

# new1 = old version of HLT_Ele30_WPTight_v + customisation
echo "Ele30 + python customisation.."
hltGetConfiguration /users/missirol/test/dev/CMSSW_15_0_0/tmp/250321_Ele30/Test03/HLT/V1 \
 --customise HLTrigger/Configuration/customizeHLTforEleWPTightRetuning.customizeHLTforEleWPTightRetuning \
 > tmp.py
edmConfigDump tmp.py > hltEle30_new_viaPython.py
rm -f tmp.py

# new2 = new version HLT_Ele30_WPTight_v implemented in ConfDB based on Laurent's config
echo "Ele30 modified via ConfDB.."
hltGetConfiguration /users/missirol/test/dev/CMSSW_15_0_0/tmp/250321_Ele30/Test03/HLT/V2 \
 > tmp.py
edmConfigDump tmp.py > hltEle30_new_viaConfDB.py
rm -f tmp.py

# show info on Paths of the GRun menu which are potentially affected by this re-tuning
echo "Info on affected Paths in GRun.. (see printInfoOnEleWPTightRetuning.txt)"
hltGetConfiguration /dev/CMSSW_15_0_0/GRun/V22 \
 --customise HLTrigger/Configuration/customizeHLTforEleWPTightRetuning.printInfoOnEleWPTightRetuning \
 > tmp.py
python3 tmp.py > printInfoOnEleWPTightRetuning.txt
rm -f tmp.py
