import FWCore.ParameterSet.Config as cms

import math

l1tHGCalTowerMapProducer = cms.EDProducer("HGCalTowerMapProducer",
    InputTriggerSums = cms.InputTag('l1tHGCalConcentratorProducer:HGCalConcentratorProcessorSelection'),
    ProcessorParameters = cms.PSet(
        ProcessorName  = cms.string('HGCalTowerMapProcessor'),
        towermap_parameters = cms.PSet(
            useLayerWeights = cms.bool(False),
            layerWeights = cms.vdouble(),
            AlgoName = cms.string('HGCalTowerMapsWrapper'),
            L1TTriggerTowerConfig = cms.PSet(
                readMappingFile = cms.bool(False),
                doNose = cms.bool(False),
                minEta = cms.double(1.479),
                maxEta = cms.double(3.0),
                minPhi = cms.double(-1*math.pi),
                maxPhi = cms.double(math.pi),
                nBinsEta = cms.int32(18),
                nBinsPhi = cms.int32(72),
                binsEta = cms.vdouble(),
                binsPhi = cms.vdouble(),
                splitModuleSum = cms.bool(False)
            )
        )
    )
)

l1tHGCalTowerMapProducerHFNose = l1tHGCalTowerMapProducer.clone(
    InputTriggerSums = 'l1tHGCalConcentratorProducerHFNose:HGCalConcentratorProcessorSelection',
    ProcessorParameters = dict(
        towermap_parameters = dict(
            L1TTriggerTowerConfig = dict(
                doNose = True,
                minEta = 3.0,
                maxEta = 4.2
            )
        )
    )
)
