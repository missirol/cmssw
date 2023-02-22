import FWCore.ParameterSet.Config as cms

import sys
import argparse

parser = argparse.ArgumentParser(prog=sys.argv[0], description='Test a simple workflow with Alpaka modules')

parser.add_argument('-a', '--accelerators', type=str, default='*',
                    help='Comma-separated string used to set process.options.accelerators (default: "*")')

parser.add_argument('-l', '--esProductLabel', type=str, default='',
                    help='Value of "appendToDataLabel" parameter of the test ESProducer (default: "")')

parser.add_argument('-s', '--useSequence', action='store_true', default=False,
                    help='Put the Alpaka EDProducer in a cms.Sequence instead of a cms.Task (default: False)')

parser.add_argument('-an', '--useAnalyzer', action='store_true', default=False,
                    help='Use EDAnalyzer to consume outputs of Alpaka EDProducer (default: False)')

parser.add_argument('-d', '--dumpPython', type=str, default=None,
                    help='Path to file containing output of process.dumpPython() (disabled by default)')

argv = sys.argv[:]
if '--' in argv: argv.remove('--')
args, unknown = parser.parse_known_args(argv)

process = cms.Process('TEST')

process.options.accelerators = args.accelerators.split(',')
print('accelerators:', process.options.accelerators.value())

if '"' in args.esProductLabel:
  args.esProductLabel = args.esProductLabel.replace('"', '')
if "'" in args.esProductLabel:
  args.esProductLabel = args.esProductLabel.replace("'", '')
print('esProductLabel:', '"'+args.esProductLabel+'"')

print('useSequence:', args.useSequence)
print('useAnalyzer:', args.useAnalyzer)

process.source = cms.Source('EmptySource')
process.maxEvents.input = 10

process.load('Configuration.StandardSequences.Accelerators_cff')
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')

process.jobConfAlpakaRcdESSource = cms.ESSource('EmptyESSource',
  recordName = cms.string('JobConfigurationAlpakaRecord'),
  iovIsRunNotTime = cms.bool(True),
  firstValid = cms.vuint32(1)
)

from RecoParticleFlow.PFRecHitProducer.pfRecHitHBHEParamsESProducer_cfi import pfRecHitHBHEParamsESProducer as _pfRecHitHBHEParamsESProducer
process.pfRecHitHBHEParamsESProducer = _pfRecHitHBHEParamsESProducer.clone(
  appendToDataLabel = args.esProductLabel
)

from RecoParticleFlow.PFRecHitProducer.testPrintPFRecHitHBHEParamsOnDevice_cfi import testPrintPFRecHitHBHEParamsOnDevice as _testPrintPFRecHitHBHEParamsOnDevice
process.testPrintPFRecHitHBHEParamsOnDevice = _testPrintPFRecHitHBHEParamsOnDevice.clone(
  pfRecHitParams = 'pfRecHitHBHEParamsESProducer:'+args.esProductLabel
)

if args.useAnalyzer:
  from RecoParticleFlow.PFRecHitProducer.testEmptyAnalyzer_cfi import testEmptyAnalyzer as _testEmptyAnalyzer
  process.testEmptyAnalyzer = _testEmptyAnalyzer.clone(
    source = 'testPrintPFRecHitHBHEParamsOnDevice'
  )

if args.useSequence:
  process.testSequence = cms.Sequence( process.testPrintPFRecHitHBHEParamsOnDevice )
  if args.useAnalyzer:
    process.testSequence += process.testEmptyAnalyzer
  process.testPath = cms.Path( process.testSequence )
else:
  process.testTask = cms.Task( process.pfRecHitHBHEParamsESProducer, process.testPrintPFRecHitHBHEParamsOnDevice )
  process.testSequence = cms.Sequence()
  if args.useAnalyzer:
    process.testSequence += process.testEmptyAnalyzer
  process.testPath = cms.Path( process.testSequence , process.testTask )

process.output = cms.OutputModule('PoolOutputModule',
  fileName = cms.untracked.string('testAlpaka.root'),
  outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_testPrintPFRecHitHBHEParamsOnDevice*_*_*',
  )
)
process.testEndPath = cms.EndPath( process.output )

# dump content of cms.Process to python file
if args.dumpPython != None:
  open(args.dumpPython, 'w').write(process.dumpPython())
