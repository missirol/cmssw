import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

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
  appendToDataLabel = '' #testProductLabel'
)

from RecoParticleFlow.PFRecHitProducer.testPrintPFRecHitHBHEParamsOnDevice_cfi import testPrintPFRecHitHBHEParamsOnDevice as _testPrintPFRecHitHBHEParamsOnDevice
process.testPrintPFRecHitHBHEParamsOnDevice = _testPrintPFRecHitHBHEParamsOnDevice.clone(
  pfRecHitParams = 'pfRecHitHBHEParamsESProducer:' #testProductLabel'
)

process.testSequence = cms.Sequence( process.testPrintPFRecHitHBHEParamsOnDevice )
process.testPath = cms.Path( process.testSequence )
