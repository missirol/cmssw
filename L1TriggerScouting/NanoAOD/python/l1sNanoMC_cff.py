import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.l1trig_cff import *
from PhysicsTools.NanoAOD.nano_cff import nanoMetadata

l1sNanoMCSequence = cms.Sequence(
    nanoMetadata
  + l1MuTable
  + l1EGTable
  + l1TauTable
  + l1JetTable
  + l1EtSumTable
)

def customiseL1ScoutNanoMCWithCaloTowers(process):

    process.l1MuTable.src = "simGmtStage2Digis"
    process.l1EGTable.src = "simCaloStage2Digis"
    process.l1TauTable.src = "simCaloStage2Digis"
    process.l1JetTable.src = "simCaloStage2Digis"
    process.l1EtSumTable.src = "simCaloStage2Digis"

    from L1TriggerScouting.Utilities.modules import L1ScoutingCaloTowerBXVecConverter
    process.l1sCaloTowerBXVec = L1ScoutingCaloTowerBXVecConverter(
        src = "simCaloStage2Layer1Digis",
    )

    from L1TriggerScouting.OnlineProcessing.modules import L1ScoutingCaloJetBXVecProducer
    process.l1sCaloJetBXVec = L1ScoutingCaloJetBXVecProducer(
        src = "l1sCaloTowerBXVec",
        akR = 0.4,
        ptMin = 5,
        applyJECs = True,
        jecFile = "L1TriggerScouting/OnlineProcessing/data/JEC_AK4CaloTowerL1S_Run3Winter25_v2.txt",
        jecPUProxyTowerMinHwEt = 1,
        jecPUProxyTowerMaxHwEt = -1,
        jecPUProxyTowerMinAbsHwEta = 0,
        jecPUProxyTowerMaxAbsHwEta = 4,
        produceSortedCaloTowers = True,
        mantissaPrecision = 10,
    )

    process.l1sCaloTowerTable = cms.EDProducer("SimpleL1ScoutingCaloTowerFlatTableProducer",
        src = cms.InputTag("l1sCaloJetBXVec:SortedCaloTowers"),
        name = cms.string("L1CaloTower"),
        minBX = cms.int32(0),
        maxBX = cms.int32(0),
        cut = cms.string("hwEt > 0"),
        doc = cms.string(""),
        extension = cms.bool(False),
        variables = cms.PSet(
            hwEt = Var("hwEt()", "int16", doc=""),
            hwEta = Var("hwEta()", "int16", doc=""),
            hwPhi = Var("hwPhi()", "int16", doc=""),
            erBits = Var("erBits()", "int16", doc=""),
            miscBits = Var("miscBits()", "int16", doc="")
        )
    )

    process.l1sCaloJetTable = cms.EDProducer("SimpleL1ScoutingCaloJetFlatTableProducer",
        src = cms.InputTag("l1sCaloJetBXVec:CaloJets"),
        name = cms.string("L1CaloJet"),
        minBX = cms.int32(0),
        maxBX = cms.int32(0),
        cut = cms.string(""),
        doc = cms.string("AK4 Jets based on CaloTowers from Calo Layer-1"),
        extension = cms.bool(False),
        variables = cms.PSet(
            pt = Var("pt()", "float", doc="jet pT", precision=10),
            eta = Var("eta()", "float", doc="jet eta", precision=10),
            phi = Var("phi()", "float", doc="jet phi", precision=10),
            mass = Var("mass()", "float", doc="jet mass", precision=10),
            energyCorr = Var("energyCorr()", "float", doc="correction factor applied to the jet-energy scale"),
            energyFracEm = Var("energyFracEm()", "float", doc="EM fraction of the jet's total energy"),
            nConst = Var("nConst()", "int", doc="number of jet constituents"),
            nConstSaturatedEnergyECAL = Var("nConstSaturatedEnergyECAL()", "uint16", doc="number of jet constituents with saturated ECAL energy"),
            nConstSaturatedEnergyHCAL = Var("nConstSaturatedEnergyHCAL()", "uint16", doc="number of jet constituents with saturated HCAL energy"),
            nConstSaturatedEnergyECALAndHCAL = Var("nConstSaturatedEnergyECALAndHCAL()", "uint16", doc="number of jet constituents with saturated energy in both ECAL and HCAL"),
        )
    )

    process.l1sNanoMCSequence += process.l1sCaloTowerBXVec
    process.l1sNanoMCSequence += process.l1sCaloJetBXVec
    process.l1sNanoMCSequence += process.l1sCaloTowerTable
    process.l1sNanoMCSequence += process.l1sCaloJetTable

    return process
