import FWCore.ParameterSet.Config as cms

def customizeHLTforCMSHLT3534(process):

    try:
        # standard SiStrips unpacker onDemand
        process.hltSiStripRawToClustersFacility.onDemand = True
        process.hltSiStripRawToClustersFacility.Clusterizer.MaxClusterSize = 16

        # SiStrips unpacker for MkFit tracks not onDemand
        process.hltSiStripRawToClustersFacilityForMkFit = process.hltSiStripRawToClustersFacility.clone(
            onDemand = False,
        )

        process.hltMeasurementTrackerEventForMkFit = process.hltMeasurementTrackerEvent.clone(
            stripClusterProducer = 'hltSiStripRawToClustersFacilityForMkFit'
        )

        # SiStrips local-reco Sequence to be used for MkFit tracks
        process.HLTDoLocalStripForMkFitSequence = cms.Sequence(
            process.hltSiStripExcludedFEDListProducer
          + process.hltSiStripRawToClustersFacilityForMkFit
          + process.hltMeasurementTrackerEventForMkFit
        )

        process.hltMeasurementTrackerEventForMkFitSerialSync = process.hltMeasurementTrackerEventSerialSync.clone(
            stripClusterProducer = 'hltSiStripRawToClustersFacilityForMkFit'
        )

        # SiStrips local-reco Sequence to be used for MkFit tracks (SerialSync variant)
        process.HLTDoLocalStripForMkFitSequenceSerialSync = cms.Sequence(
            process.hltSiStripExcludedFEDListProducer
          + process.hltSiStripRawToClustersFacilityForMkFit
          + process.hltMeasurementTrackerEventForMkFitSerialSync
        )
    except: pass

    # update SiStrips clusters used by the relevant modules producing or consuming MkFit tracks
    try:
        process.hltDoubletRecoveryClustersRefRemoval.stripClusters = 'hltSiStripRawToClustersFacilityForMkFit'
        process.hltDoubletRecoveryMaskedMeasurementTrackerEvent.src = 'hltMeasurementTrackerEventForMkFit'
        process.hltDoubletRecoveryPFlowPixelClusterCheck.ClusterCollectionLabel = 'hltMeasurementTrackerEventForMkFit'
    except: pass

    try:
        process.hltDoubletRecoveryClustersRefRemovalSerialSync.stripClusters = 'hltSiStripRawToClustersFacilityForMkFit'
        process.hltDoubletRecoveryMaskedMeasurementTrackerEventSerialSync.src = 'hltMeasurementTrackerEventForMkFitSerialSync'
        process.hltDoubletRecoveryPFlowPixelClusterCheckSerialSync.ClusterCollectionLabel = 'hltMeasurementTrackerEventForMkFitSerialSync'
    except: pass

    try:
        process.hltSiStripRecHits.ClusterProducer = 'hltSiStripRawToClustersFacilityForMkFit'
        process.hltIter0PFlowCkfTrackCandidatesMkFitSiStripHits.clusters = 'hltSiStripRawToClustersFacilityForMkFit'
        process.hltIter0PFlowCtfWithMaterialTracks.MeasurementTrackerEvent = 'hltMeasurementTrackerEventForMkFit'
    except: pass

    try:
        process.hltIter0PFlowCtfWithMaterialTracksSerialSync.MeasurementTrackerEvent = 'hltMeasurementTrackerEventForMkFitSerialSync'
    except: pass

    try:
        process.hltDisplacedhltIter4ClustersRefRemovalForTau.stripClusters = 'hltSiStripRawToClustersFacilityForMkFit'
        process.hltDisplacedhltIter4MaskedMeasurementTrackerEventForTau.src = 'hltMeasurementTrackerEventForMkFit'
        process.hltDisplacedhltIter4PFlowPixelLessClusterCheckForTau.ClusterCollectionLabel = 'hltMeasurementTrackerEventForMkFit'
    except: pass

    # The Sequence named HLTTrackReconstructionForPFNoMu contains the Sequences used for
    #  - (1) pixel local reconstruction plus pixel tracking,
    #  - (2) strips local reconstruction for Iter0+Iter2 tracks, and
    #  - (3) Iter0+Iter2 tracks, used for the standard ParticleFlow Sequence.
    #
    # Below, the following changes are made.
    #  - (a) HLTDoLocalStripSequence is replaced by HLTDoLocalStripForMkFitSequence
    #        in the Sequences which include HLTIterativeTrackingIter02 (MkFit tracks).
    #  - (b) A SerialSync variant of the Sequence HLTTrackReconstructionForPFNoMu is introduced,
    #        i.e. HLTTrackReconstructionForPFNoMuSerialSync.
    #  - (c) Other Sequences and Paths which contain explicitly the Sequences corresponding to (1)+(2)+(3)
    #        are updated in order to use HLTTrackReconstructionForPFNoMu(SerialSync) instead.
    #        This way, for future updates similar to this one, it will be sufficient
    #        to update the content of the HLTTrackReconstructionForPFNoMu(SerialSync) Sequences,
    #        as opposed to touching many other Sequences and Paths (like in this case).
    try:
        process.HLTTrackReconstructionForPFNoMu = cms.Sequence(
            process.HLTDoLocalPixelSequence
          + process.HLTRecopixelvertexingSequence
          + process.HLTDoLocalStripForMkFitSequence
          + process.HLTIterativeTrackingIter02
        )

        process.HLTTrackReconstructionForPFNoMuSerialSync = cms.Sequence(
            process.HLTDoLocalPixelSequenceSerialSync
          + process.HLTRecopixelvertexingSequenceSerialSync
          + process.HLTDoLocalStripForMkFitSequenceSerialSync
          + process.HLTIterativeTrackingIter02SerialSync
        )
    except: pass

    # Redefine the relevant Sequences and Paths so they use HLTTrackReconstructionForPFNoMu(SerialSync).
    try:
        # this Sequence is used by AlCa_PFJet40_v (and many other Paths)
        process.HLTTrackReconstructionForPF = cms.Sequence(
            process.HLTTrackReconstructionForPFNoMu
          + process.hltPFMuonMerging
          + process.hltMuonLinks
          + process.hltMuons
        )
    except: pass

    try:
        # this Sequence is used by AlCa_PFJet40_CPUOnly_v
        process.HLTTrackReconstructionForPFSerialSync = cms.Sequence(
            process.HLTTrackReconstructionForPFNoMuSerialSync
          + process.hltPFMuonMergingSerialSync
          + process.hltMuonLinksSerialSync
          + process.hltMuonsSerialSync
        )
    except: pass

    try:
        # this Sequence is used by HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24_v (and many other Paths)
        process.HLTIterativeTrackingIter04ForTau = cms.Sequence(
            process.HLTTrackReconstructionForPFNoMu
          + process.HLTIterativeTrackingIteration4ForTau
          + process.hltIter4MergedWithIter0ForTau
        )

        process.HLTTrackReconstructionForPFDispl = cms.Sequence(
            process.HLTIterativeTrackingIter04ForTau
          + process.hltPFMuonMergingForDisplTau
          + process.hltMuonLinksForDisplTau
          + process.hltMuonsForDisplTau
        )
    except: pass

    try:
        # this Sequence is used by HLT_DoubleMu3_TkMu_DsTau3Mu_v (and many other Paths)
        process.HLTTrackerMuonSequenceLowPt = cms.Sequence(
            process.HLTTrackReconstructionForPFNoMu
          + process.HLTL3muonrecoNocandSequence
          + process.hltDiMuonMergingIter01TkMu
          + process.hltDiMuonLinksIter01TkMuMerge
          + process.hltGlbTrkMuonsLowPtIter01Merge
          + process.hltGlbTrkMuonLowPtIter01MergeCands
        )
    except: pass

    try:
        # this Sequence is used by HLT_ZeroBias_Beamspot_v (and many other Paths)
        process.HLTTrackingForBeamSpot = cms.Sequence(
            process.HLTPreAK4PFJetsRecoSequence
          + process.HLTL2muonrecoSequence
          + process.HLTL3muonrecoSequence
          + process.HLTTrackReconstructionForPFNoMu
          + process.hltPFMuonMerging
        )
    except: pass

    # replace HLTDoLocalStripSequence with HLTDoLocalStripForMkFitSequence
    # in the Paths that used the former to produce MkFit tracks (HLTIterativeTrackingIter02)
    for pathName, path in process.paths_().items():
        if not pathName.startswith('AlCa_IsoTrackHBHE_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sHTTMultiJet
          + process.hltPreAlCaIsoTrackHBHE
          + process.HLTTrackReconstructionForPFNoMu
          + process.hltMergedTracksSelector
          + process.hltMergedTracksSelectorFilter
          + process.HLTDoFullUnpackingEgammaEcalSequence
          + process.HLTDoLocalHcalSequence
          + process.HLTL2muonrecoSequence
          + process.HLTL3muonrecoSequence
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_MET105_IsoTrk50_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sETM90ToETM150
          + process.hltPreMET105IsoTrk50
          + process.HLTRecoMETSequence
          + process.hltMET105
          + process.HLTRecoJetSequenceAK4PrePF
          + process.HLTTrackReconstructionForPFNoMu
          + process.hltDeDxEstimatorProducer
          + process.hltTrk50Filter
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_MET120_IsoTrk50_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sETM80ToETM150
          + process.hltPreMET120IsoTrk50
          + process.HLTRecoMETSequence
          + process.hltMET120
          + process.HLTRecoJetSequenceAK4PrePF
          + process.HLTTrackReconstructionForPFNoMu
          + process.hltDeDxEstimatorProducer
          + process.hltTrk50Filter
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_IsoTrk200_L1SingleMuShower_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sMuShowerOneNominal
          + process.hltPreIsoTrk200L1SingleMuShower
          + process.HLTTrackReconstructionForPFNoMu
          + process.hltDeDxEstimatorProducer
          + process.hltTrk200MuonEndcapFilter
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_IsoTrk400_L1SingleMuShower_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sMuShowerOneNominal
          + process.hltPreIsoTrk400L1SingleMuShower
          + process.HLTTrackReconstructionForPFNoMu
          + process.hltDeDxEstimatorProducer
          + process.hltTrk400MuonEndcapFilter
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_PFMET105_IsoTrk50_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sETM90ToETM150
          + process.hltPrePFMET105IsoTrk50
          + process.HLTRecoMETSequence
          + process.hltMET75
          + process.HLTRecoJetSequenceAK4PrePF
          + process.HLTTrackReconstructionForPFNoMu
          + process.hltDeDxEstimatorProducer
          + process.hltTrk50Filter
          + process.HLTAK4PFJetsSequence
          + process.hltPFMETProducer
          + process.hltPFMET105
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_IsoMu24_HLTTracking_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sSingleMu22
          + process.hltPreIsoMu24HLTTracking
          + process.hltL1fL1sMu22L1Filtered0
          + process.HLTL2muonrecoSequence
          + cms.ignore(process.hltL2fL1sSingleMu22L1f0L2Filtered10Q)
          + process.HLTL3muonrecoSequence
          + cms.ignore(process.hltL1fForIterL3L1fL1sMu22L1Filtered0)
          + process.hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q
          + process.HLTMu24IsolationSequence
          + process.hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered
          + process.HLTTrackReconstructionForPFNoMu
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_HLTTracking_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sDoubleMu125to157
          + process.hltPreMu17TrkIsoVVLMu8TrkIsoVVLDZHLTTracking
          + process.hltL1fL1sDoubleMu155L1Filtered0
          + process.HLTL2muonrecoSequence
          + cms.ignore(process.hltL2pfL1sDoubleMu155L1f0L2PreFiltered0)
          + cms.ignore(process.hltL2fL1sDoubleMu155L1f0L2Filtered10OneMu)
          + process.HLTL3muonrecoSequence
          + cms.ignore(process.hltL1fForIterL3L1fL1sDoubleMu155L1Filtered0)
          + process.hltL3fL1DoubleMu155fPreFiltered8
          + process.hltL3fL1DoubleMu155fFiltered17
          + process.HLTL3muontrkisovvlSequence
          + process.hltDiMuon178RelTrkIsoVVLFiltered
          + process.hltDiMuon178RelTrkIsoVVLFilteredDzFiltered0p2
          + process.HLTTrackReconstructionForPFNoMu
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_HLTTracking_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltL1sDoubleMu125to157
          + process.hltPreMu17TrkIsoVVLMu8TrkIsoVVLDZMass3p8HLTTracking
          + process.hltL1fL1sDoubleMu155L1Filtered0
          + process.HLTL2muonrecoSequence
          + cms.ignore(process.hltL2pfL1sDoubleMu155L1f0L2PreFiltered0)
          + cms.ignore(process.hltL2fL1sDoubleMu155L1f0L2Filtered10OneMu)
          + process.HLTL3muonrecoSequence
          + cms.ignore(process.hltL1fForIterL3L1fL1sDoubleMu155L1Filtered0)
          + process.hltL3fL1DoubleMu155fPreFiltered8
          + process.hltL3fL1DoubleMu155fFiltered17
          + process.HLTL3muontrkisovvlSequence
          + process.hltDiMuon178RelTrkIsoVVLFiltered
          + process.hltDiMuon178RelTrkIsoVVLFilteredDzFiltered0p2
          + process.hltDiMuon178Mass3p8Filtered
          + process.HLTTrackReconstructionForPFNoMu
          + process.HLTEndSequence
        ))

    for pathName, path in process.paths_().items():
        if not pathName.startswith('MC_ReducedIterativeTracking_v'): continue
        setattr(process, pathName, cms.Path(
            process.HLTBeginSequence
          + process.hltPreMCReducedIterativeTracking
          + process.HLTTrackReconstructionForPFNoMu
          + process.HLTEndSequence
        ))

    # add hltSiStripRawToClustersFacilityForMkFit to the EventContent of
    # all the OutputModules in which hltSiStripRawToClustersFacility is kept
    #  - Streams: DQM, HLTMonitor, CosmicHLTMonitor.
    for outMod in process.outputModules_().values():
        extra_keeps = []
        for outCmd in outMod.outputCommands:
            if 'hltSiStripRawToClustersFacility' in outCmd:
                extra_keeps += [outCmd.replace('hltSiStripRawToClustersFacility', 'hltSiStripRawToClustersFacilityForMkFit')]
        outMod.outputCommands += extra_keeps

    return process
