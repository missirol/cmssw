#!/usr/bin/env python3
import sys

if len(sys.argv) != 2:
  raise Exception('specify 1 cmd-line argument: the type of HLT menu (must be "GRun", "HIon", "PIon", or "PRef")')

hltMenuType = sys.argv[1]

if hltMenuType == 'GRun':
    from HLTrigger.Configuration.HLT_GRun_cff import cms,fragment

elif hltMenuType == 'HIon':
    from HLTrigger.Configuration.HLT_HIon_cff import cms,fragment

elif hltMenuType == 'PIon':
    from HLTrigger.Configuration.HLT_PIon_cff import cms,fragment

elif hltMenuType == 'PRef':
    from HLTrigger.Configuration.HLT_PRef_cff import cms,fragment

else:
  raise Exception('invalid type of HLT menu (must be "GRun", "HIon", "PIon", or "PRef")')

for dsetName in fragment.datasets.parameterNames_():
  print('Dataset_'+dsetName)

#for ttt in GRun HIon PIon PRef; do
#  echo "" >> "${ttt}".txt
#  ./printDatasets.py "${ttt}" >> "${ttt}".txt
#  echo "" >> online_"${ttt,,}".txt
#  ./printDatasets.py "${ttt}" >> online_"${ttt,,}".txt
#done
#unset ttt
