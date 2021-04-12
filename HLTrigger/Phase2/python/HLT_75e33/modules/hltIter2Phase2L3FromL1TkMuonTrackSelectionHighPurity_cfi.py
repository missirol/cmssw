import FWCore.ParameterSet.Config as cms

hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter2Phase2L3FromL1TkMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter2Phase2L3FromL1TkMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks")
)
