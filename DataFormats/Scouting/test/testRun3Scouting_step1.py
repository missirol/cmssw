import FWCore.ParameterSet.Config as cms

## CLI parser
import argparse
import sys

parser = argparse.ArgumentParser(
  prog = 'cmsRun '+sys.argv[0]+' --',
  description = 'Configuration file to test I/O of Scouting collections.',
  formatter_class = argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('-t', '--nThreads', type = int, help = 'Number of threads',
                    default = 1)

parser.add_argument('-s', '--nStreams', type = int, help = 'Number of EDM streams',
                    default = 0)

parser.add_argument('-i', '--inputFiles', nargs = '+', help = 'List of EDM input files',
                    default = ['/store/data/Run2022A/ScoutingPFRun3/RAW/v1/000/352/565/00000/d95fdb9d-14fb-471f-8abf-a3bf31faa804.root'])

parser.add_argument('-n', '--maxEvents', type = int, help = 'Number of input events',
                    default = 100)

parser.add_argument('-o', '--outputFile', type = str, help = 'Path to output EDM file in ROOT format',
                    default = 'testDataFormatsScouting_output.root')

parser.add_argument('--wantSummary', action = 'store_true', help = 'Value of process.options.wantSummary',
                    default = False)

argv = sys.argv[:]
if '--' in argv:
    argv.remove('--')
args, unknown = parser.parse_known_args(argv)

## Process
process = cms.Process('TEST')

process.options.numberOfThreads = args.nThreads
process.options.numberOfStreams = args.nStreams
process.options.wantSummary = args.wantSummary

process.maxEvents.input = args.maxEvents

# Source (EDM input)
process.source = cms.Source('PoolSource',
  fileNames = cms.untracked.vstring(args.inputFiles),
  inputCommands = cms.untracked.vstring(
    'drop *',
    'keep *Scouting*_*_*_*',
  )
)

# MessageLogger (Service)
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Output module
process.testOutput = cms.OutputModule('PoolOutputModule',
  fileName = cms.untracked.string( args.outputFile )
)

# EndPath
process.testEndPath = cms.EndPath( process.testOutput )
