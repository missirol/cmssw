#!/bin/bash

hltLabel=testCMSHLT3534

hltConfigFromDB --configName /dev/CMSSW_15_0_0/HLT > tmp.py

edmConfigDump tmp.py > "${hltLabel}"_hlt0.py

cat <<@EOF >> tmp.py

from HLTrigger.Configuration.customizeHLTforCMSHLT3534 import customizeHLTforCMSHLT3534
process = customizeHLTforCMSHLT3534(process)
@EOF

edmConfigDump tmp.py > "${hltLabel}"_hlt1.py

rm -f tmp.py
unset hltLabel
