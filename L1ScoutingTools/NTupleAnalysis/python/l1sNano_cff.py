import FWCore.ParameterSet.Config as cms

from L1TriggerScouting.OnlineProcessing.l1sAK4CTJets0Emu_cfi import l1sAK4CTJets0Emu
from L1TriggerScouting.OnlineProcessing.l1sAK4CTJets1Emu_cfi import l1sAK4CTJets1Emu

from L1ScoutingTools.NanoAOD.l1sNanoTables_cff import *

from PhysicsTools.NanoAOD.l1trig_cff import l1TablesTask
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles

from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import prunedGenParticles

l1sNanoTask = cms.Task(nanoMetadata)

l1sNanoSequence = cms.Sequence(l1sNanoTask)

l1EmulExtraObjsTask = cms.Task(
    l1sAK4CTJets0Emu,
    l1sAK4CTJets1Emu,
)

def customiseNanoForL1ScoutCaloTowersMC(process):
    process.l1sNanoTask.add(process.l1EmulExtraObjsTask)
    process.l1sNanoTask.add(process.l1EmulCaloLayer1NanoTask)
    process.l1sNanoTask.add(process.l1EmulObjTablesTask)
    return process
