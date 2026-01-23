import FWCore.ParameterSet.Config as cms

from L1ScoutingTools.NanoAOD.l1sNanoTables_cff import *

from L1ScoutingTools.Reconstruction.l1sAK4CTJets0Emu_cfi import l1sAK4CTJets0Emu
from L1ScoutingTools.Reconstruction.l1sAK4CTJets1Emu_cfi import l1sAK4CTJets1Emu

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
    process = customizeNanoForL1ScoutGEN(process)
    return process

def customizeNanoForL1ScoutGEN(process):
    # delete unnecessary GEN-related products
    for foo in [
        'ak4GenJetsChargedOnly',
        'genFilterTable',
        'genIso',
        'genParticles2HepMC',
        'genParticles2HepMCHiggsVtx',
        'genParticlesForJetsCharged',
        'HTXSCategoryTable',
        'mergedGenParticles',
        'particleLevel',
        'rivetProducerHTXS',
        'tautagger',
        'trackGenJetAK4Table',
    ]:
        if hasattr(process, foo):
            delattr(process, foo)

    try:
        # remove sequences related to particle-level information
        process.nanogenSequence.remove(process.particleLevelSequence)
        process.nanogenSequence.remove(process.particleLevelTablesSequence)
    except:
        pass

    try:
        # restrict list of GenParticles in the output
        process.prunedGenParticles.src = 'genParticles'
        process.genParticleTable.src = 'finalGenParticles'
        process.genParticleTablesTask.add(process.prunedGenParticles)
        process.genParticleTablesTask.add(process.finalGenParticles)
    except:
        pass

    try:
        # add pileup-related information
        process.puTable.src = 'addPileupInfo'
        process.puTable.savePUDensityVars = False
        process.puTable.pvsrc = ''
        process.puTable.zbins = []
        if not process.nanogenSequence.contains(process.puTable):
            process.nanogenSequence.insert(0, process.puTable)
    except:
        pass

    return process
