cmsRun "${CMSSW_BASE}"/src/HLTrigger/HLTfilters/test/triggerResultsFilter_testOperatorMasking.py &> tmp.log && grep -m15 'TrigReport ' tmp.log
