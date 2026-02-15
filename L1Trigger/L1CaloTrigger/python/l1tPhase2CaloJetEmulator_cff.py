from L1Trigger.L1CaloTrigger.l1tPhase2CaloJetEmulator_cfi import *

from L1Trigger.L1THGCal.L1TTriggerTowerConfig_energySplit_cfi import L1TTriggerTowerConfig_energySplit as _L1TTriggerTowerConfig_energySplit
from L1Trigger.L1THGCal.l1tHGCalTowerMapProducer_cfi import *
from L1Trigger.L1THGCal.l1tHGCalTowerProducer_cfi import *

# Add HGCal tower producers for energy split towers
# Based on custom_towers_energySplit in L1Trigger/L1THGCal/python/customTowers.py
l1tHGCalEnergySplitTowerMapProducer = l1tHGCalTowerMapProducer.clone(
    ProcessorParameters = dict(
        towermap_parameters = dict(
            L1TTriggerTowerConfig = _L1TTriggerTowerConfig_energySplit.clone()
        )
    )
)

l1tHGCalEnergySplitTowerProducer = l1tHGCalTowerProducer.clone(
    InputTowerMaps = ("l1tHGCalEnergySplitTowerMapProducer","HGCalTowerMapProcessor")
)

l1tHGCalEnergySplitTowersTask = cms.Task(
    l1tHGCalEnergySplitTowerMapProducer,
    l1tHGCalEnergySplitTowerProducer
)

# Use energy split towers in calo jet/tau emulator
l1tPhase2CaloJetEmulator.hgcalTowers = ("l1tHGCalEnergySplitTowerProducer","HGCalTowerProcessor")

l1tCaloJetsTausTask = cms.Task(
    l1tHGCalEnergySplitTowersTask,
    l1tPhase2CaloJetEmulator
)
