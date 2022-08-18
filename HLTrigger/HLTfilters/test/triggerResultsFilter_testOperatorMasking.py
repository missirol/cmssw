import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')
process.options.wantSummary = True
process.options.numberOfThreads = 1

process.load('FWCore.MessageService.MessageLogger_cfi')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    reportEvery = cms.untracked.int32(1), # every!
#    limit = cms.untracked.int32(-1)       # no limit!
#)
process.MessageLogger.cerr.FwkReport.reportEvery = 1 # only report every 100th event start
#process.MessageLogger.cerr_stats.threshold = 'INFO' # also info in statistics

# define the Prescaler service, and set the scale factors
process.PrescaleService = cms.Service('PrescaleService',
  lvl1DefaultLabel = cms.string('any'),
  lvl1Labels = cms.vstring('any'),
  prescaleTable = cms.VPSet(
    cms.PSet(
      pathName  = cms.string('HLT_OnlyEven'),
      prescales = cms.vuint32( 2 )
    ),
    cms.PSet(
      pathName  = cms.string('HLT_OnlyOdd'),
      prescales = cms.vuint32( 2 )
    )
  )
)

# define an empty source, and ask for 100 events
process.source = cms.Source('EmptySource')
process.maxEvents.input = 3

process.hltPreOnlyEven = cms.EDFilter('HLTPrescaler', offset = cms.uint32(0))
process.hltPreOnlyOdd  = cms.EDFilter('HLTPrescaler', offset = cms.uint32(1))
process.hltFilterAlwaysTrue  = cms.EDFilter('HLTBool', result = cms.bool(True))
process.hltFilterAlwaysFalse = cms.EDFilter('HLTBool', result = cms.bool(False))

process.HLT_OnlyEven = cms.Path(process.hltPreOnlyEven)
process.HLT_OnlyOdd = cms.Path(process.hltPreOnlyOdd)
process.HLT_AlwaysTrue = cms.Path(process.hltFilterAlwaysTrue)
process.HLT_AlwaysFalse = cms.Path(process.hltFilterAlwaysFalse)

# define the TriggerResultsFilters based on the status of the previous paths
from HLTrigger.HLTfilters.triggerResultsFilter_cfi import triggerResultsFilter as _trigResFilter
_triggerResultsFilter = _trigResFilter.clone( usePathStatus = True, l1tResults = '', throw = True )

def _addPath(process, name : str = '', triggerConditions = []):
  filterName = 'hltFilterCheck'+name
  setattr(process, filterName, _triggerResultsFilter.clone( triggerConditions = triggerConditions ))
  setattr(process, 'Check_'+name, cms.Path( getattr(process, filterName) ))
  print('TrigReport | Check_'+name+' = ', getattr(process, filterName).triggerConditions.value())
  return process

#process = _addPath(process, '1', ['(HLT_Only* AND HLT_AlwaysFalse) MASKING HLT_AlwaysFalse'])
#process = _addPath(process, '2', ['(HLT_Only* AND HLT_AlwaysFalse) MASKING HLT_Always*'])
#process = _addPath(process, '3', ['(HLT_OnlyEven OR HLT_OnlyOdd) MASKING (HLT_OnlyEven OR HLT_OnlyOdd)'])
#process = _addPath(process, '4', ['(HLT_OnlyEven OR HLT_OnlyOdd) MASKING (HLT_OnlyEven AND HLT_OnlyOdd)'])
#process = _addPath(process, '5', ['(HLT_OnlyEven OR HLT_OnlyOdd) MASKING HLT_OnlyEven MASKING HLT_OnlyOdd'])
#process = _addPath(process, '6', ['(HLT_OnlyEven OR HLT_OnlyOdd) MASKING HLT_Only* MASKING HLT_OnlyOdd'])
process = _addPath(process, '7', ['(HLT_OnlyEven OR HLT_OnlyOdd) MASKING (HLT_OnlyEven MASKING HLT_OnlyOdd)'])

# define an EndPath to analyze all other path results
process.hltTrigReport = cms.EDAnalyzer( 'HLTrigReport',
  HLTriggerResults = cms.InputTag( 'TriggerResults', '', '@currentProcess' )
)
process.HLTAnalyzerEndpath = cms.EndPath( process.hltTrigReport )
