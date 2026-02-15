import FWCore.ParameterSet.Config as cms

RPCConeBuilder = cms.ESProducer("RPCConeBuilder",
    towerBeg = cms.int32(0),
    towerEnd = cms.int32(16)
)

rpcconesrc = cms.ESSource("EmptyESSource",
    recordName = cms.string('L1RPCConeBuilderRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

# remove RPCConeBuilder-related ES modules when the modifier stage2L1Trigger is enabled
def _removeRPCConeBuilderESModules(process):
    if hasattr(process, 'RPCConeBuilder'):
        del process.RPCConeBuilder
    if hasattr(process, 'rpcconesrc'):
        del process.rpcconesrc

from Configuration.Eras.Modifier_stage2L1Trigger_cff import stage2L1Trigger
removeRPCConeBuilderESModules_ = stage2L1Trigger.makeProcessModifier( _removeRPCConeBuilderESModules )
