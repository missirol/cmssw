#!/usr/bin/env python3

# hltConfigFromDB --configName /dev/CMSSW_12_4_0/HIon > tsg.py
from tsg import cms,process as tsg
# hltConfigFromDB --configName /users/chenyi/PbPb2022/IonTest/HLT/V17 > hion.py
from hion import process as hin

import os

tsgPaths = sorted(list((foo[:foo.rfind('_v')]+'_v\*' if '_v' in foo else foo) for foo in tsg.paths_() if not foo.startswith('Dataset_')))
hinPaths = sorted(list((foo[:foo.rfind('_v')]+'_v\*' if '_v' in foo else foo) for foo in hin.paths_() if not foo.startswith('Dataset_')))

tsgFinalPaths = sorted(list(tsg.finalpaths_()))
hinFinalPaths = sorted(list(hin.finalpaths_()))

#hinPaths = [foo.replace('HIDataRawPrime','HITestRawPrime').replace('HIData','HITestRaw') for foo in hinPaths]
#hinFinalPaths = [foo.replace('HIDataRawPrime','HITestRawPrime').replace('HIData','HITestRaw') for foo in hinFinalPaths]

for [tsg_paths,hin_paths] in [
  [tsgPaths, hinPaths],
  [tsgFinalPaths, hinFinalPaths],
]:
  for hin_p in hin_paths:
    if hin_p not in tsg_paths:
      raise RuntimeError('missing: '+hin_p)

  for tsg_p in tsg_paths:
    if tsg_p not in hin_paths:
#      os.system('sed -i "s|'+tsg_p+'||g" HLTrigger/Configuration/tables/HIon.txt')
      os.system('sed -i "s|'+tsg_p+'|#'+tsg_p+'|g" HLTrigger/Configuration/tables/online_hion.txt')
