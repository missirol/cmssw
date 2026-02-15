import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCal.l1tHGCalTowerMapProducer_cfi import l1tHGCalTowerMapProducer as _l1tHGCalTowerMapProducer, l1tHGCalTowerMapProducerHFNose as _l1tHGCalTowerMapProducerHFNose

l1tHGCalTowerProducer = cms.EDProducer("HGCalTowerProducer",
    InputTowerMaps = cms.InputTag('l1tHGCalTowerMapProducer:HGCalTowerMapProcessor'), 
    InputTriggerCells = cms.InputTag('l1tHGCalBackEndLayer1Producer:HGCalBackendLayer1Processor2DClustering'),
    ProcessorParameters = cms.PSet(
        ProcessorName = cms.string('HGCalTowerProcessor'),
        includeTrigCells = cms.bool(False),
        towermap_parameters = _l1tHGCalTowerMapProducer.ProcessorParameters.towermap_parameters.clone()
    )
)

l1tHGCalTowerProducerHFNose = l1tHGCalTowerProducer.clone(
    InputTowerMaps = 'l1tHGCalTowerMapProducerHFNose:HGCalTowerMapProcessor',
    InputTriggerCells = 'l1tHGCalBackEndLayer1ProducerHFNose:HGCalBackendLayer1Processor2DClustering',
    ProcessorParameters = dict(
        towermap_parameters = _l1tHGCalTowerMapProducerHFNose.ProcessorParameters.towermap_parameters.clone()
    )
)

