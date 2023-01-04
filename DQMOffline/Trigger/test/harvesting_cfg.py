import FWCore.ParameterSet.Config as cms

# VarParsing
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')
options.register('nThreads', 4, options.multiplicity.singleton, options.varType.int, 'number of threads')
options.register('nStreams', 0, options.multiplicity.singleton, options.varType.int, 'number of streams')
options.setDefault('inputFiles', ['file:DQMIO.root'])
options.parseArguments()

# Process
process = cms.Process('HARVESTING')

process.options.numberOfThreads = options.nThreads
process.options.numberOfStreams = options.nStreams
process.options.numberOfConcurrentLuminosityBlocks = 1

# Source (DQM input)
process.source = cms.Source('DQMRootSource',
  fileNames = cms.untracked.vstring(options.inputFiles)
)

# DQMStore (Service)
process.load('DQMServices.Core.DQMStore_cfi')

# MessageLogger (Service)
process.load('FWCore.MessageLogger.MessageLogger_cfi')

# Output module (file in ROOT format)
from DQMServices.Components.DQMFileSaver_cfi import dqmSaver as _dqmSaver
process.dqmSaver = _dqmSaver.clone(
  workflow = '/DQMOffline/Trigger/'+process.name_()
)

# EndPath
process.endp = cms.EndPath(
  process.dqmSaver
)
