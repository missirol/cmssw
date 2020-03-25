import FWCore.ParameterSet.Config as cms

from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *

pfClusterRefsForJetsHCAL = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterHCAL'),
    particleType = cms.string('pi+')
)

pfClusterRefsForJetsECAL = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterECAL'),
    particleType = cms.string('pi+')
)

pfClusterRefsForJetsHF = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterHF'),
    particleType = cms.string('pi+')
)

pfClusterRefsForJetsHO = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterHO'),
    particleType = cms.string('pi+')
)

pfClusterRefsForJets = cms.EDProducer("PFClusterRefCandidateMerger",
    src = cms.VInputTag("pfClusterRefsForJetsHCAL", "pfClusterRefsForJetsECAL", "pfClusterRefsForJetsHF", "pfClusterRefsForJetsHO")
)

pfClusterRefsForJets_stepTask = cms.Task(
   particleFlowClusterTask,
   pfClusterRefsForJetsHCAL,
   pfClusterRefsForJetsECAL,
   pfClusterRefsForJetsHF,
   pfClusterRefsForJetsHO,
   pfClusterRefsForJets
)
pfClusterRefsForJets_step = cms.Sequence(pfClusterRefsForJets_stepTask)

### Phase-2, HGCal

pfClusterRefsForJetsHGCal = cms.EDProducer('PFClusterRefCandidateProducer',
  src = cms.InputTag('particleFlowClusterHGCal'),
  particleType = cms.string('pi+')
)

_pfClusterRefsForJets_phase2_hgcal_src = pfClusterRefsForJets.src + ['pfClusterRefsForJetsHGCal']

from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
phase2_hgcal.toModify( pfClusterRefsForJets, src = _pfClusterRefsForJets_phase2_hgcal_src )

_pfClusterRefsForJets_stepTask_phase2_hgcal = pfClusterRefsForJets_stepTask.copy()
_pfClusterRefsForJets_stepTask_phase2_hgcal.add(pfClusterRefsForJetsHGCal)

phase2_hgcal.toReplaceWith( pfClusterRefsForJets_stepTask, _pfClusterRefsForJets_stepTask_phase2_hgcal )
