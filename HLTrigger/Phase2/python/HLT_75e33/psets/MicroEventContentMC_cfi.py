import FWCore.ParameterSet.Config as cms

MicroEventContentMC = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_slimmedPhotons_*_*', 
        'keep *_slimmedOOTPhotons_*_*', 
        'keep *_slimmedElectrons_*_*', 
        'keep *_slimmedMuons_*_*', 
        'keep *_slimmedTaus_*_*', 
        'keep *_slimmedTausBoosted_*_*', 
        'keep *_slimmedCaloJets_*_*', 
        'keep *_slimmedJets_*_*', 
        'keep recoBaseTagInfosOwned_slimmedJets_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'drop recoBaseTagInfosOwned_slimmedJetsAK8_*_*', 
        'keep *_slimmedJetsPuppi_*_*', 
        'keep *_slimmedMETs_*_*', 
        'keep *_slimmedMETsNoHF_*_*', 
        'keep *_slimmedMETsPuppi_*_*', 
        'keep *_slimmedSecondaryVertices_*_*', 
        'keep *_slimmedLambdaVertices_*_*', 
        'keep *_slimmedKshortVertices_*_*', 
        'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_SubJets_*', 
        'keep recoPhotonCores_reducedEgamma_*_*', 
        'keep recoGsfElectronCores_reducedEgamma_*_*', 
        'keep recoConversions_reducedEgamma_*_*', 
        'keep recoSuperClusters_reducedEgamma_*_*', 
        'keep recoCaloClusters_reducedEgamma_*_*', 
        'keep EcalRecHitsSorted_reducedEgamma_*_*', 
        'keep recoGsfTracks_reducedEgamma_*_*', 
        'keep HBHERecHitsSorted_reducedEgamma_*_*', 
        'drop *_*_caloTowers_*', 
        'drop *_*_pfCandidates_*', 
        'drop *_*_genJets_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlineSlimmedPrimaryVertices_*_*', 
        'keep patPackedCandidates_packedPFCandidates_*_*', 
        'keep *_isolatedTracks_*_*', 
        'keep *_oniaPhotonCandidates_*_*', 
        'keep *_bunchSpacingProducer_*_*', 
        'keep double_fixedGridRhoAll__*', 
        'keep double_fixedGridRhoFastjetAll__*', 
        'keep double_fixedGridRhoFastjetAllTmp__*', 
        'keep double_fixedGridRhoFastjetAllCalo__*', 
        'keep double_fixedGridRhoFastjetCentral_*_*', 
        'keep double_fixedGridRhoFastjetCentralCalo__*', 
        'keep double_fixedGridRhoFastjetCentralChargedPileUp__*', 
        'keep double_fixedGridRhoFastjetCentralNeutral__*', 
        'keep *_slimmedPatTrigger_*_*', 
        'keep patPackedTriggerPrescales_patTrigger__*', 
        'keep patPackedTriggerPrescales_patTrigger_l1max_*', 
        'keep patPackedTriggerPrescales_patTrigger_l1min_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_gtStage2Digis__*', 
        'keep *_gmtStage2Digis_Muon_*', 
        'keep *_caloStage2Digis_Jet_*', 
        'keep *_caloStage2Digis_Tau_*', 
        'keep *_caloStage2Digis_EGamma_*', 
        'keep *_caloStage2Digis_EtSum_*', 
        'keep *_TriggerResults_*_HLT', 
        'keep *_TriggerResults_*_*', 
        'keep patPackedCandidates_lostTracks_*_*', 
        'keep HcalNoiseSummary_hcalnoise__*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*', 
        'keep recoForwardProtons_ctppsProtons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep *_prefiringweight_*_*', 
        'keep *_packedPFCandidates_hcalDepthEnergyFractions_*', 
        'drop *_packedPFCandidates_hcalDepthEnergyFractions_*', 
        'keep patPackedGenParticles_packedGenParticles_*_*', 
        'keep recoGenParticles_prunedGenParticles_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_*_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep recoGenParticles_genPUProtons_*_*', 
        'keep *_slimmedGenJetsFlavourInfos_*_*', 
        'keep *_slimmedGenJets__*', 
        'keep *_slimmedGenJetsAK8__*', 
        'keep *_slimmedGenJetsAK8SoftDropSubJets__*', 
        'keep *_genMetTrue_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep GenRunInfoProduct_*_*_*', 
        'keep *_genParticles_xyz0_*', 
        'keep *_genParticles_t0_*', 
        'keep PileupSummaryInfos_slimmedAddPileupInfo_*_*', 
        'keep L1GtTriggerMenuLite_l1GtTriggerMenuLite__*', 
        'keep *_slimmedElectronsFromMultiCl_*_*', 
        'keep *_slimmedPhotonsFromMultiCl_*_*', 
        'keep *_offlineSlimmedPrimaryVertices4D_*_*'
    )
)