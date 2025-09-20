import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

process.options.numberOfThreads = 1
process.options.numberOfStreams = 0
process.options.wantSummary = False

process.maxEvents.input = -1

# MessageLogger
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Input source
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        '/store/data/Run2025F/ZeroBias/RAW/v1/000/397/209/00000/7477b158-4f18-4e7e-ac0c-a8f603e9fe65.root'
    )
)

# EventData modules
process.eventNumberFilter = cms.EDFilter("HLTL1NumberFilter",
    rawInput = cms.InputTag( "rawDataCollector" ),
    period = cms.uint32( 0 ),
    invert = cms.bool( False ),
    fedId = cms.int32( 0 ),
    useTCDSEventNumber = cms.bool( False )
)

# Path definition
process.Path1 = cms.Path( process.eventNumberFilter )
