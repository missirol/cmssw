import FWCore.ParameterSet.Config as cms

softPFMuonComputer = cms.ESProducer("MuonTaggerESProducer",
    gbrForestLabel = cms.string('btag_SoftPFMuon_BDT'),
    ipSign = cms.string('any'),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFMuon_BDT.weights.xml.gz')
)
