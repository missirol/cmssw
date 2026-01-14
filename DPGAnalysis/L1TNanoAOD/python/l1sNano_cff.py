import FWCore.ParameterSet.Config as cms

from DPGAnalysis.L1TNanoAOD.l1sNanoTables_cff import *
from PhysicsTools.NanoAOD.l1trig_cff import l1TablesTask
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles
from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import prunedGenParticles

l1sNanoTask = cms.Task(nanoMetadata)

l1sNanoSequence = cms.Sequence(l1sNanoTask)

def customiseNanoForL1ScoutCaloTowersMC(process):
    process.l1sNanoTask.add(process.l1EmulCaloLayer1NanoTask)
    process.l1sNanoTask.add(process.l1EmulObjTablesTask)

    # delete unnecessary GEN-related products (if present)
    for genLabel in [
        'trackGenJetAK4Table',
        'HTXSCategoryTable',
        'genParticlesForJetsCharged',
        'ak4GenJetsChargedOnly',
        'genParticles2HepMC',
        'genParticles2HepMCHiggsVtx',
        'genParticlesForJetsCharged',
        'ak4GenJetsChargedOnly',
        'tautagger',
        'rivetProducerHTXS',
        'rivetLeptonTable',
        'rivetMetTable',
        'rivetPhotonTable',
    ]:
        if hasattr(process, genLabel):
            delattr(process, genLabel)

    # restrict list of GenParticles in the output
    try:
       process.genParticleTablesTask.add(process.prunedGenParticles)
       process.genParticleTablesTask.add(process.finalGenParticles)
       process.prunedGenParticles.src = 'genParticles'
       process.genParticleTable.src = 'finalGenParticles'
    except:
       pass

    return process
