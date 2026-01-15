import FWCore.ParameterSet.Config as cms

from L1TriggerScouting.OnlineProcessing.L1TCaloTowerAKJetProducer import L1TCaloTowerAKJetProducer

l1sAK4CTJetsEmu = L1TCaloTowerAKJetProducer(
    src = 'simCaloStage2Layer1Digis',
    bxMin = -2,
    bxMax = 2,
    towerMinHwPt = 1,
    rParam = 0.4,
    jetPtMin = 5
)
