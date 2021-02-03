import FWCore.ParameterSet.Config as cms

REPACKRAWSIMEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'drop *', 
        'drop FEDRawDataCollection_*_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep FEDRawDataCollection_rawDataReducedFormat_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep *_hltScoutingTrackPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_simMuonCSCDigis_*_*', 
        'keep *_simMuonRPCDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonGEMDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonME0Digis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep *_simHcalUnsuppressedDigis_*_*', 
        'keep *_mix_EETimeDigi_*', 
        'keep *_mix_EBTimeDigi_*', 
        'keep *_simEcalUnsuppressedDigis_*_*', 
        'keep *_simHGCalUnsuppressedDigis_EE_*', 
        'keep *_simHGCalUnsuppressedDigis_HEfront_*', 
        'keep *_simHGCalUnsuppressedDigis_HEback_*', 
        'keep *_mix_MergedCaloTruth_*', 
        'keep *_mix_FTLBarrel_*', 
        'keep *_mix_FTLEndcap_*', 
        'keep *_mix_InitialVertices_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_ak*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'drop FEDRawDataCollection_source_*_*', 
        'drop FEDRawDataCollection_rawDataCollector_*_*'
    ),
    splitLevel = cms.untracked.int32(0)
)