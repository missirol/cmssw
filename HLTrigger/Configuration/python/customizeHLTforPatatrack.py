import FWCore.ParameterSet.Config as cms

from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA
from HLTrigger.Configuration.common import *
from Configuration.Eras.Modifier_run3_common_cff import run3_common


# force the SwitchProducerCUDA choice to pick a specific backend: True for offloading to a gpu, False for running on cpu
def forceGpuOffload(status = True):
    import HeterogeneousCore.CUDACore.SwitchProducerCUDA
    HeterogeneousCore.CUDACore.SwitchProducerCUDA._cuda_enabled_cached = bool(status)


# reset the SwitchProducerCUDA choice to pick a backend depending on the availability of a supported gpu
def resetGpuOffload():
    import HeterogeneousCore.CUDACore.SwitchProducerCUDA
    HeterogeneousCore.CUDACore.SwitchProducerCUDA._cuda_enabled_cached = None
    HeterogeneousCore.CUDACore.SwitchProducerCUDA._switch_cuda()


# customisation for running the Patatrack reconstruction, common parts
def customiseCommon(process):

    # Services

    if 'CUDAService' not in process.__dict__:
        process.load("HeterogeneousCore.CUDAServices.CUDAService_cfi")

    if 'MessageLogger' in process.__dict__:
        process.MessageLogger.CUDAService = cms.untracked.PSet()

#    # NVProfilerService is broken in CMSSW 12.0.x and later
#    if 'NVProfilerService' in process.__dict__:
#        process.load("HeterogeneousCore.CUDAServices.NVProfilerService_cfi")


    # Paths and EndPaths

    # the hltGetConditions module would force gpu-specific ESProducers to run even if no supported gpu is present
    if 'hltGetConditions' in process.__dict__:
        del process.hltGetConditions

    # produce a boolean to track if the events ar being processed on gpu (true) or cpu (false)
    if 'statusOnGPU' not in process.__dict__:
        process.statusOnGPU = SwitchProducerCUDA(
            cpu  = cms.EDProducer("BooleanProducer", value = cms.bool(False)),
            cuda = cms.EDProducer("BooleanProducer", value = cms.bool(True))
        )

    if 'statusOnGPUFilter' not in process.__dict__:
        process.statusOnGPUFilter = cms.EDFilter("BooleanFilter",
            src = cms.InputTag("statusOnGPU")
        )

    if 'Status_OnCPU' in process.__dict__:
        replace_with(process.Status_OnCPU, cms.Path(process.statusOnGPU + ~process.statusOnGPUFilter))
    else:
        process.Status_OnCPU = cms.Path(process.statusOnGPU + ~process.statusOnGPUFilter)
        if process.schedule is not None:
            process.schedule.append(process.Status_OnCPU)

    if 'Status_OnGPU' in process.__dict__:
        replace_with(process.Status_OnGPU, cms.Path(process.statusOnGPU + process.statusOnGPUFilter))
    else:
        process.Status_OnGPU = cms.Path(process.statusOnGPU + process.statusOnGPUFilter)
        if process.schedule is not None:
            process.schedule.append(process.Status_OnGPU)

    # make the ScoutingCaloMuonOutput endpath compatible with using Tasks in the Scouting paths
    if 'hltOutputScoutingCaloMuon' in process.__dict__ and not 'hltPreScoutingCaloMuonOutputSmart' in process.__dict__:
        process.hltPreScoutingCaloMuonOutputSmart = cms.EDFilter( "TriggerResultsFilter",
            l1tIgnoreMaskAndPrescale = cms.bool( False ),
            l1tResults = cms.InputTag( "" ),
            hltResults = cms.InputTag( 'TriggerResults','','@currentProcess' ),
            triggerConditions = process.hltOutputScoutingCaloMuon.SelectEvents.SelectEvents,
            throw = cms.bool( True )
        )
        insert_modules_after(process, process.hltPreScoutingCaloMuonOutput, process.hltPreScoutingCaloMuonOutputSmart)

    # make the ScoutingPFOutput endpath compatible with using Tasks in the Scouting paths
    if 'hltOutputScoutingPF' in process.__dict__ and not 'hltPreScoutingPFOutputSmart' in process.__dict__:
        process.hltPreScoutingPFOutputSmart = cms.EDFilter( "TriggerResultsFilter",
            l1tIgnoreMaskAndPrescale = cms.bool( False ),
            l1tResults = cms.InputTag( "" ),
            hltResults = cms.InputTag( 'TriggerResults','','@currentProcess' ),
            triggerConditions = process.hltOutputScoutingPF.SelectEvents.SelectEvents,
            throw = cms.bool( True )
        )
        insert_modules_after(process, process.hltPreScoutingPFOutput, process.hltPreScoutingPFOutputSmart)


    # done
    return process


# customisation for running the "Patatrack" pixel local reconstruction
def customisePixelLocalReconstruction(process):

    if not 'HLTDoLocalPixelSequence' in process.__dict__:
        return process


    # FIXME replace the Sequences with empty ones to avoid expanding them during the (re)definition of Modules and EDAliases
    process.HLTDoLocalPixelSequence = cms.Sequence()


    # Event Setup
    if 'siPixelGainCalibrationForHLTGPU' not in process.__dict__:
        process.load("CalibTracker.SiPixelESProducers.siPixelGainCalibrationForHLTGPU_cfi")                 # this should be used only on GPUs, will crash otherwise

    if 'siPixelROCsStatusAndMappingWrapperESProducer' not in process.__dict__:
        process.load("CalibTracker.SiPixelESProducers.siPixelROCsStatusAndMappingWrapperESProducer_cfi")    # this should be used only on GPUs, will crash otherwise

    if 'PixelCPEFastESProducer' not in process.__dict__:
        process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEFastESProducer_cfi")


    # Modules and EDAliases
    # referenced in HLTDoLocalPixelTask

    if 'hltOnlineBeamSpotToCUDA' not in process.__dict__:
        # transfer the beamspot to the gpu
        from RecoVertex.BeamSpotProducer.offlineBeamSpotToCUDA_cfi import offlineBeamSpotToCUDA as _offlineBeamSpotToCUDA
        process.hltOnlineBeamSpotToCUDA = _offlineBeamSpotToCUDA.clone(
            src = "hltOnlineBeamSpot"
        )

    if 'hltSiPixelClustersCUDA' not in process.__dict__:
        # reconstruct the pixel digis and clusters on the gpu
        _hltSiPixelClustersLegacyLabel = 'hltSiPixelClustersLegacy' if 'hltSiPixelClustersLegacy' in process.__dict__ else 'hltSiPixelClusters'

        from RecoLocalTracker.SiPixelClusterizer.siPixelRawToClusterCUDA_cfi import siPixelRawToClusterCUDA as _siPixelRawToClusterCUDA
        process.hltSiPixelClustersCUDA = _siPixelRawToClusterCUDA.clone(
            # use the same thresholds as the legacy module
            clusterThreshold_layer1 = getattr(process, _hltSiPixelClustersLegacyLabel).ClusterThreshold_L1,
            clusterThreshold_otherLayers = getattr(process, _hltSiPixelClustersLegacyLabel).ClusterThreshold,
        )
        # use the pixel channel calibrations scheme for Run 3
        run3_common.toModify(process.hltSiPixelClustersCUDA, isRun2 = False)

    # copy the pixel digis errors to the host
    if 'hltSiPixelDigiErrorsSoA' not in process.__dict__:
        from EventFilter.SiPixelRawToDigi.siPixelDigiErrorsSoAFromCUDA_cfi import siPixelDigiErrorsSoAFromCUDA as _siPixelDigiErrorsSoAFromCUDA
        process.hltSiPixelDigiErrorsSoA = _siPixelDigiErrorsSoAFromCUDA.clone(
            src = "hltSiPixelClustersCUDA"
        )

    # copy the pixel digis (except errors) and clusters to the host
    if 'hltSiPixelDigisSoA' not in process.__dict__:
        from EventFilter.SiPixelRawToDigi.siPixelDigisSoAFromCUDA_cfi import siPixelDigisSoAFromCUDA as _siPixelDigisSoAFromCUDA
        process.hltSiPixelDigisSoA = _siPixelDigisSoAFromCUDA.clone(
            src = "hltSiPixelClustersCUDA"
        )


    # SwitchProducer: hltSiPixelDigis
    if not isinstance(process.hltSiPixelDigis, SwitchProducerCUDA):

        if 'hltSiPixelDigisLegacy' not in process.__dict__:
            # reconstruct the pixel digis on the cpu
            process.hltSiPixelDigisLegacy = process.hltSiPixelDigis.clone()

        # SwitchProducer wrapping a subset of the legacy pixel digis producer, or the conversion of the pixel digis errors to the legacy format
        process.hltSiPixelDigis = SwitchProducerCUDA(
            # legacy producer
            cpu = cms.EDAlias(
                hltSiPixelDigisLegacy = cms.VPSet(
                    cms.PSet(type = cms.string("DetIdedmEDCollection")),
                    cms.PSet(type = cms.string("SiPixelRawDataErroredmDetSetVector")),
                    cms.PSet(type = cms.string("PixelFEDChanneledmNewDetSetVector"))
                )
            )
        )

    elif not hasattr(process.hltSiPixelDigis, 'cpu'):
        raise Exception('unsupported configuration: "process.hltSiPixelDigis" is a SwitchProducerCUDA, but does not have a "cpu" branch')

    if not hasattr(process.hltSiPixelDigis, 'cuda'):
        # conversion from SoA to legacy format
        from EventFilter.SiPixelRawToDigi.siPixelDigiErrorsFromSoA_cfi import siPixelDigiErrorsFromSoA as _siPixelDigiErrorsFromSoA
        process.hltSiPixelDigis.cuda = _siPixelDigiErrorsFromSoA.clone(
            digiErrorSoASrc = "hltSiPixelDigiErrorsSoA",
            UsePhase1 = True
        )


    # SwitchProducer: hltSiPixelClusters
    if not isinstance(process.hltSiPixelClusters, SwitchProducerCUDA):

        if 'hltSiPixelClustersLegacy' not in process.__dict__:
            # reconstruct the pixel clusters on the cpu
            process.hltSiPixelClustersLegacy = process.hltSiPixelClusters.clone(
                src = "hltSiPixelDigisLegacy"
            )

        # SwitchProducer wrapping a subset of the legacy pixel cluster producer, or the conversion of the pixel digis (except errors) and clusters to the legacy format
        process.hltSiPixelClusters = SwitchProducerCUDA(
            # legacy producer
            cpu = cms.EDAlias(
                hltSiPixelClustersLegacy = cms.VPSet(
                    cms.PSet(type = cms.string("SiPixelClusteredmNewDetSetVector"))
                )
            )
        )

    elif not hasattr(process.hltSiPixelClusters, 'cpu'):
        raise Exception('unsupported configuration: "process.hltSiPixelClusters" is a SwitchProducerCUDA, but does not have a "cpu" branch')

    if not hasattr(process.hltSiPixelClusters, 'cuda'):
        # conversion from SoA to legacy format
        from RecoLocalTracker.SiPixelClusterizer.siPixelDigisClustersFromSoA_cfi import siPixelDigisClustersFromSoA as _siPixelDigisClustersFromSoA
        process.hltSiPixelClusters.cuda = _siPixelDigisClustersFromSoA.clone(
            src = "hltSiPixelDigisSoA",
            produceDigis = False,
            storeDigis = False,
            # use the same thresholds as the legacy module
            clusterThreshold_layer1 = process.hltSiPixelClustersLegacy.ClusterThreshold_L1,
            clusterThreshold_otherLayers = process.hltSiPixelClustersLegacy.ClusterThreshold
        )


    # reconstruct the pixel rechits on the gpu
    if 'hltSiPixelRecHitsCUDA' not in process.__dict__:
        from RecoLocalTracker.SiPixelRecHits.siPixelRecHitCUDA_cfi import siPixelRecHitCUDA as _siPixelRecHitCUDA
        process.hltSiPixelRecHitsCUDA = _siPixelRecHitCUDA.clone(
            src = "hltSiPixelClustersCUDA",
            beamSpot = "hltOnlineBeamSpotToCUDA"
        )

    # cpu only: produce the pixel rechits in SoA and legacy format, from the legacy clusters
    if 'hltSiPixelRecHitSoA' not in process.__dict__:
        from RecoLocalTracker.SiPixelRecHits.siPixelRecHitSoAFromLegacy_cfi import siPixelRecHitSoAFromLegacy as _siPixelRecHitSoAFromLegacy
        process.hltSiPixelRecHitSoA = _siPixelRecHitSoAFromLegacy.clone(
            src = "hltSiPixelClusters",
            beamSpot = "hltOnlineBeamSpot",
            convertToLegacy = True
        )


    # SwitchProducer: hltSiPixelRecHits
    if not isinstance(process.hltSiPixelRecHits, SwitchProducerCUDA):
        # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA
        process.hltSiPixelRecHits = SwitchProducerCUDA(
            # legacy producer
            cpu = cms.EDAlias(
                hltSiPixelRecHitSoA = cms.VPSet(
                    cms.PSet(type = cms.string("SiPixelRecHitedmNewDetSetVector")),
                    cms.PSet(type = cms.string("uintAsHostProduct"))
                )
            )
        )

    elif not hasattr(process.hltSiPixelRecHits, 'cpu'):
        raise Exception('unsupported configuration: "process.hltSiPixelRecHits" is a SwitchProducerCUDA, but does not have a "cpu" branch')

    if not hasattr(process.hltSiPixelRecHits, 'cuda'):
        # conversion from SoA to legacy format
        from RecoLocalTracker.SiPixelRecHits.siPixelRecHitFromCUDA_cfi import siPixelRecHitFromCUDA as _siPixelRecHitFromCUDA
        process.hltSiPixelRecHits.cuda = _siPixelRecHitFromCUDA.clone(
            pixelRecHitSrc = "hltSiPixelRecHitsCUDA",
            src = "hltSiPixelClusters",
        )


    # Tasks and Sequences
    if 'HLTDoLocalPixelTask' not in process.__dict__:
        process.HLTDoLocalPixelTask = cms.Task(
            process.hltOnlineBeamSpotToCUDA,   # transfer the beamspot to the gpu
            process.hltSiPixelClustersCUDA,    # reconstruct the pixel digis and clusters on the gpu
            process.hltSiPixelRecHitsCUDA,     # reconstruct the pixel rechits on the gpu
            process.hltSiPixelDigisSoA,        # copy the pixel digis (except errors) and clusters to the host
            process.hltSiPixelDigiErrorsSoA,   # copy the pixel digis errors to the host
            process.hltSiPixelDigisLegacy,     # legacy pixel digis producer
            process.hltSiPixelDigis,           # SwitchProducer wrapping a subset of the legacy pixel digis producer, or the conversion of the pixel digis errors from SoA
            process.hltSiPixelClustersLegacy,  # legacy pixel cluster producer
            process.hltSiPixelClusters,        # SwitchProducer wrapping a subset of the legacy pixel cluster producer, or the conversion of the pixel digis (except errors) and clusters from SoA
            process.hltSiPixelClustersCache,   # legacy module, used by the legacy pixel quadruplet producer
            process.hltSiPixelRecHitSoA,       # pixel rechits on cpu, in SoA & legacy format
            process.hltSiPixelRecHits,         # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA
        )

    elif not isinstance(process.HLTDoLocalPixelTask, cms.Task):
        raise Exception('unsupported configuration: "process.HLTDoLocalPixelTask" already exists, but it is not a Task')

    # redefine HLTDoLocalPixelSequence (it was emptied at the start of this function)
    process.HLTDoLocalPixelSequence = cms.Sequence(process.HLTDoLocalPixelTask)


    # workaround for old version of AlCa_LumiPixelsCounts paths
    for AlCaPathName in ['AlCa_LumiPixelsCounts_Random_v1', 'AlCa_LumiPixelsCounts_ZeroBias_v1']:
        if AlCaPathName in process.__dict__:
            AlCaPath = getattr(process, AlCaPathName)
            # replace hltSiPixelDigis+hltSiPixelClusters with HLTDoLocalPixelSequence
            hasSiPixelDigis, hasSiPixelClusters = False, False
            for (itemLabel, itemName) in AlCaPath.directDependencies():
                if itemLabel != 'modules': continue
                if itemName == 'hltSiPixelDigis': hasSiPixelDigis = True
                elif itemName == 'hltSiPixelClusters': hasSiPixelClusters = True
            if hasSiPixelDigis and hasSiPixelClusters:
                AlCaPath.remove(process.hltSiPixelClusters)
                AlCaPath.replace(process.hltSiPixelDigis, process.HLTDoLocalPixelSequence)


    # done
    return process


# customisation for running the "Patatrack" pixel track reconstruction
def customisePixelTrackReconstruction(process):

    if not 'HLTRecoPixelTracksSequence' in process.__dict__:
        return process

    hasHLTPixelVertexReco = 'HLTRecopixelvertexingSequence' in process.__dict__

    # FIXME replace the Sequences with empty ones to avoid expanding them during the (re)definition of Modules and EDAliases
    process.HLTRecoPixelTracksSequence = cms.Sequence()
    if hasHLTPixelVertexReco:
        process.HLTRecopixelvertexingSequence = cms.Sequence()


    # Modules and EDAliases
    # referenced in process.HLTRecoPixelTracksTask

    # SwitchProducer: hltPixelTracksSoA
    if not ('hltPixelTracksSoA' in process.__dict__ and isinstance(process.hltPixelTracksSoA, SwitchProducerCUDA)):
        # build pixel ntuplets and pixel tracks in SoA format on gpu
        from RecoPixelVertexing.PixelTriplets.pixelTracksCUDA_cfi import pixelTracksCUDA as _pixelTracksCUDA
        process.hltPixelTracksCUDA = _pixelTracksCUDA.clone(
            pixelRecHitSrc = "hltSiPixelRecHitsCUDA",
            onGPU = True,
            idealConditions = False,
        )
        # use quality cuts tuned for Run-2 ideal conditions for all Run-3 workflows
        run3_common.toModify(process.hltPixelTracksCUDA, idealConditions = True)

        process.hltPixelTracksSoA = SwitchProducerCUDA(
            # build pixel ntuplets and pixel tracks in SoA format on cpu
            cpu = process.hltPixelTracksCUDA.clone(
                pixelRecHitSrc = "hltSiPixelRecHitSoA",
                onGPU = False,
            )
        )

    elif hasattr(process.hltPixelTracksSoA, 'cpu'):
        # if cpu branch of SwitchProducerCUDA exists, take hltPixelTracksCUDA (gpu)
        # from hltPixelTracksSoA.cpu (cpu) to enforce same configuration parameters
        process.hltPixelTracksCUDA = process.hltPixelTracksSoA.cpu.clone(
            pixelRecHitSrc = "hltSiPixelRecHitsCUDA",
            onGPU = True
        )

    else:
        raise Exception('unsupported configuration: "process.hltPixelTracksSoA" is a SwitchProducerCUDA, but does not have a "cpu" branch')

    # if cuda branch of SwitchProducerCUDA does not exist, add it
    if not hasattr(process.hltPixelTracksSoA, 'cuda'):
        # transfer the pixel tracks in SoA format to cpu
        from RecoPixelVertexing.PixelTrackFitting.pixelTracksSoA_cfi import pixelTracksSoA as _pixelTracksSoA
        process.hltPixelTracksSoA.cuda = _pixelTracksSoA.clone(
            src = "hltPixelTracksCUDA"
        )


    # convert the pixel tracks from SoA to legacy format
    if process.hltPixelTracks.type_() != 'PixelTrackProducerFromSoA':
        from RecoPixelVertexing.PixelTrackFitting.pixelTrackProducerFromSoA_cfi import pixelTrackProducerFromSoA as _pixelTrackProducerFromSoA
        process.hltPixelTracks = _pixelTrackProducerFromSoA.clone(
            beamSpot = "hltOnlineBeamSpot",
            pixelRecHitLegacySrc = "hltSiPixelRecHits",
            trackSrc = "hltPixelTracksSoA",
        )


    # referenced in process.HLTRecopixelvertexingTask
    if hasHLTPixelVertexReco:

        # SwitchProducer: hltPixelVerticesSoA
        if not ('hltPixelVerticesSoA' in process.__dict__ and isinstance(process.hltPixelVerticesSoA, SwitchProducerCUDA)):
            # build pixel vertices in SoA format on gpu
            from RecoPixelVertexing.PixelVertexFinding.pixelVerticesCUDA_cfi import pixelVerticesCUDA as _pixelVerticesCUDA
            process.hltPixelVerticesCUDA = _pixelVerticesCUDA.clone(
                pixelTrackSrc = "hltPixelTracksCUDA",
                onGPU = True
            )

            # build or transfer pixel vertices in SoA format on cpu
            process.hltPixelVerticesSoA = SwitchProducerCUDA(
                # build pixel vertices in SoA format on cpu
                cpu = _pixelVerticesCUDA.clone(
                    pixelTrackSrc = "hltPixelTracksSoA",
                    onGPU = False
                )
            )

        elif hasattr(process.hltPixelVerticesSoA, 'cpu'):
            # if cpu branch of SwitchProducerCUDA exists, take hltPixelVerticesCUDA (gpu)
            # from hltPixelVerticesSoA.cpu (cpu) to enforce same configuration parameters
            process.hltPixelVerticesCUDA = process.hltPixelVerticesSoA.cpu.clone(
                pixelTrackSrc = "hltPixelTracksCUDA",
                onGPU = True
            )

        else:
            raise Exception('unsupported configuration: "process.hltPixelVerticesSoA" is a SwitchProducerCUDA, but does not have a "cpu" branch')

        # if cuda branch of SwitchProducerCUDA does not exist, add it
        if not hasattr(process.hltPixelVerticesSoA, 'cuda'):
            # transfer the pixel vertices in SoA format to cpu
            from RecoPixelVertexing.PixelVertexFinding.pixelVerticesSoA_cfi import pixelVerticesSoA as _pixelVerticesSoA
            process.hltPixelVerticesSoA.cuda = _pixelVerticesSoA.clone(
                src = "hltPixelVerticesCUDA"
            )


        # convert the pixel vertices from SoA to legacy format
        if process.hltPixelTracks.type_() != 'PixelVertexProducerFromSoA':
            from RecoPixelVertexing.PixelVertexFinding.pixelVertexFromSoA_cfi import pixelVertexFromSoA as _pixelVertexFromSoA
            process.hltPixelVertices = _pixelVertexFromSoA.clone(
                src = "hltPixelVerticesSoA",
                TrackCollection = "hltPixelTracks",
                beamSpot = "hltOnlineBeamSpot",
            )


    # Tasks and Sequences

    if 'HLTRecoPixelTracksTask' not in process.__dict__:
        process.HLTRecoPixelTracksTask = cms.Task(
            process.hltPixelTracksTrackingRegions,         # from the original sequence
            process.hltPixelTracksCUDA,                    # pixel ntuplets on gpu, in SoA format
            process.hltPixelTracksSoA,                     # pixel ntuplets on cpu, in SoA format
            process.hltPixelTracks,                        # pixel tracks on cpu, in legacy format
        )

    elif not isinstance(process.HLTRecoPixelTracksTask, cms.Task):
        raise Exception('unsupported configuration: "process.HLTRecoPixelTracksTask" already exists, but it is not a Task')

    # redefine HLTRecoPixelTracksSequence (it was emptied at the start of this function)
    process.HLTRecoPixelTracksSequence = cms.Sequence(process.HLTRecoPixelTracksTask)

    if hasHLTPixelVertexReco:

        if 'HLTRecopixelvertexingTask' not in process.__dict__:
            process.HLTRecopixelvertexingTask = cms.Task(
                process.HLTRecoPixelTracksTask,
                process.hltPixelVerticesCUDA,              # pixel vertices on gpu, in SoA format
                process.hltPixelVerticesSoA,               # pixel vertices on cpu, in SoA format
                process.hltPixelVertices,                  # pixel vertices on cpu, in legacy format
                process.hltTrimmedPixelVertices,           # from the original sequence
            )

        elif not isinstance(process.HLTRecopixelvertexingTask, cms.Task):
            raise Exception('unsupported configuration: "process.HLTRecopixelvertexingTask" already exists, but it is not a Task')

        # redefine HLTRecopixelvertexingSequence (it was emptied at the start of this function)
        process.HLTRecopixelvertexingSequence = cms.Sequence(
            process.hltPixelTracksFitter +             # not used here, kept for compatibility with legacy sequences
            process.hltPixelTracksFilter,              # not used here, kept for compatibility with legacy sequences
            process.HLTRecopixelvertexingTask,
        )


    # done
    return process


# customisation for offloading the ECAL local reconstruction via CUDA if a supported gpu is present
def customiseEcalLocalReconstruction(process):

    hasHLTEcalPreshowerSeq = any(seq in process.__dict__ for seq in ['HLTDoFullUnpackingEgammaEcalMFSequence', 'HLTDoFullUnpackingEgammaEcalSequence'])
    if not (hasHLTEcalPreshowerSeq or 'HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence' in process.__dict__):
        return process

    # FIXME replace the Sequences with empty ones to avoid expanding them during the (re)definition of Modules and EDAliases
    process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence = cms.Sequence()
    if hasHLTEcalPreshowerSeq:
        process.HLTDoFullUnpackingEgammaEcalMFSequence = cms.Sequence()
        process.HLTDoFullUnpackingEgammaEcalSequence = cms.Sequence()


    # Event Setup

    if 'ecalElectronicsMappingGPUESProducer' not in process.__dict__:
        process.load("EventFilter.EcalRawToDigi.ecalElectronicsMappingGPUESProducer_cfi")

    if 'ecalGainRatiosGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalGainRatiosGPUESProducer_cfi")

    if 'ecalPedestalsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalPedestalsGPUESProducer_cfi")

    if 'ecalPulseCovariancesGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalPulseCovariancesGPUESProducer_cfi")

    if 'ecalPulseShapesGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalPulseShapesGPUESProducer_cfi")

    if 'ecalSamplesCorrelationGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalSamplesCorrelationGPUESProducer_cfi")

    if 'ecalTimeBiasCorrectionsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalTimeBiasCorrectionsGPUESProducer_cfi")

    if 'ecalTimeCalibConstantsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalTimeCalibConstantsGPUESProducer_cfi")

    if 'ecalMultifitParametersGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalMultifitParametersGPUESProducer_cfi")

    if 'ecalRechitADCToGeVConstantGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalRechitADCToGeVConstantGPUESProducer_cfi")

    if 'ecalRechitChannelStatusGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalRechitChannelStatusGPUESProducer_cfi")

    if 'ecalIntercalibConstantsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalIntercalibConstantsGPUESProducer_cfi")

    if 'ecalLaserAPDPNRatiosGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalLaserAPDPNRatiosGPUESProducer_cfi")

    if 'ecalLaserAPDPNRatiosRefGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalLaserAPDPNRatiosRefGPUESProducer_cfi")

    if 'ecalLaserAlphasGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalLaserAlphasGPUESProducer_cfi")

    if 'ecalLinearCorrectionsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalLinearCorrectionsGPUESProducer_cfi")

    if 'ecalRecHitParametersGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.EcalRecProducers.ecalRecHitParametersGPUESProducer_cfi")


    # Modules and EDAliases

    if 'hltEcalDigisGPU' not in process.__dict__:
        # ECAL unpacker running on gpu
        from EventFilter.EcalRawToDigi.ecalRawToDigiGPU_cfi import ecalRawToDigiGPU as _ecalRawToDigiGPU
        process.hltEcalDigisGPU = _ecalRawToDigiGPU.clone()


    # SwitchProducer: hltEcalDigis
    if not isinstance(process.hltEcalDigis, SwitchProducerCUDA):

        if 'hltEcalDigisLegacy' not in process.__dict__:
            process.hltEcalDigisLegacy = process.hltEcalDigis.clone()
        else:
            raise Exception('unsupported configuration: "process.hltEcalDigis" is not a SwitchProducerCUDA, but "process.hltEcalDigisLegacy" already exists')

        # SwitchProducer wrapping the legacy ECAL unpacker or the ECAL digi converter from SoA format on gpu to legacy format on cpu
        process.hltEcalDigis = SwitchProducerCUDA(
            # legacy producer
            cpu = cms.EDAlias(
                hltEcalDigisLegacy = cms.VPSet(
                    cms.PSet(type = cms.string("EBDigiCollection")),
                    cms.PSet(type = cms.string("EEDigiCollection")),
                    cms.PSet(type = cms.string("EBDetIdedmEDCollection")),
                    cms.PSet(type = cms.string("EEDetIdedmEDCollection")),
                    cms.PSet(type = cms.string("EBSrFlagsSorted")),
                    cms.PSet(type = cms.string("EESrFlagsSorted")),
                    cms.PSet(type = cms.string("EcalElectronicsIdedmEDCollection"), fromProductInstance = cms.string("EcalIntegrityBlockSizeErrors")),
                    cms.PSet(type = cms.string("EcalElectronicsIdedmEDCollection"), fromProductInstance = cms.string("EcalIntegrityTTIdErrors")),
                    cms.PSet(type = cms.string("EcalElectronicsIdedmEDCollection"), fromProductInstance = cms.string("EcalIntegrityZSXtalIdErrors")),
                    cms.PSet(type = cms.string("EcalPnDiodeDigisSorted")),
                    cms.PSet(type = cms.string("EcalPseudoStripInputDigisSorted"), fromProductInstance = cms.string("EcalPseudoStripInputs")),
                    cms.PSet(type = cms.string("EcalTriggerPrimitiveDigisSorted"), fromProductInstance = cms.string("EcalTriggerPrimitives")),
                )
            )
        )

    if not hasattr(process.hltEcalDigis, 'cuda'):
        # convert ECAL digis from SoA format on gpu to legacy format on cpu
        from EventFilter.EcalRawToDigi.ecalCPUDigisProducer_cfi import ecalCPUDigisProducer as _ecalCPUDigisProducer
        process.hltEcalDigis.cuda = _ecalCPUDigisProducer.clone(
            digisInLabelEB = ("hltEcalDigisGPU", "ebDigis"),
            digisInLabelEE = ("hltEcalDigisGPU", "eeDigis"),
            produceDummyIntegrityCollections = True
        )

    if 'hltEcalUncalibRecHitGPU' not in process.__dict__:
        # ECAL multifit running on gpu
        from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitProducerGPU_cfi import ecalUncalibRecHitProducerGPU as _ecalUncalibRecHitProducerGPU
        process.hltEcalUncalibRecHitGPU = _ecalUncalibRecHitProducerGPU.clone(
            digisLabelEB = ("hltEcalDigisGPU", "ebDigis"),
            digisLabelEE = ("hltEcalDigisGPU", "eeDigis"),
            shouldRunTimingComputation = False
        )

    if 'hltEcalUncalibRecHitSoA' not in process.__dict__:
        # copy the ECAL uncalibrated rechits from gpu to cpu in SoA format
        from RecoLocalCalo.EcalRecProducers.ecalCPUUncalibRecHitProducer_cfi import ecalCPUUncalibRecHitProducer as _ecalCPUUncalibRecHitProducer
        process.hltEcalUncalibRecHitSoA = _ecalCPUUncalibRecHitProducer.clone(
            recHitsInLabelEB = ("hltEcalUncalibRecHitGPU", "EcalUncalibRecHitsEB"),
            recHitsInLabelEE = ("hltEcalUncalibRecHitGPU", "EcalUncalibRecHitsEE"),
        )

    # SwitchProducer: hltEcalUncalibRecHit
    if not isinstance(process.hltEcalUncalibRecHit, SwitchProducerCUDA):
        # SwitchProducer wrapping the legacy ECAL uncalibrated rechits producer or a converter from SoA to legacy format
        process.hltEcalUncalibRecHit = SwitchProducerCUDA(
            # legacy producer
            cpu = process.hltEcalUncalibRecHit
        )

    if not hasattr(process.hltEcalUncalibRecHit, 'cuda'):
        # convert the ECAL uncalibrated rechits from SoA to legacy format
        from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitConvertGPU2CPUFormat_cfi import ecalUncalibRecHitConvertGPU2CPUFormat as _ecalUncalibRecHitConvertGPU2CPUFormat
        process.hltEcalUncalibRecHit.cuda = _ecalUncalibRecHitConvertGPU2CPUFormat.clone(
            recHitsLabelGPUEB = ("hltEcalUncalibRecHitSoA", "EcalUncalibRecHitsEB"),
            recHitsLabelGPUEE = ("hltEcalUncalibRecHitSoA", "EcalUncalibRecHitsEE"),
        )

    # Reconstructing the ECAL calibrated rechits on gpu works, but is extremely slow.
    # Disable it for the time being, until the performance has been addressed.
    """
    from RecoLocalCalo.EcalRecProducers.ecalRecHitGPU_cfi import ecalRecHitGPU as _ecalRecHitGPU
    process.hltEcalRecHitGPU = _ecalRecHitGPU.clone(
        uncalibrecHitsInLabelEB = ("hltEcalUncalibRecHitGPU","EcalUncalibRecHitsEB"),
        uncalibrecHitsInLabelEE = ("hltEcalUncalibRecHitGPU","EcalUncalibRecHitsEE"),
    )

    from RecoLocalCalo.EcalRecProducers.ecalCPURecHitProducer_cfi import ecalCPURecHitProducer as _ecalCPURecHitProducer
    process.hltEcalRecHitSoA = _ecalCPURecHitProducer.clone(
        recHitsInLabelEB = ("hltEcalRecHitGPU", "EcalRecHitsEB"),
        recHitsInLabelEE = ("hltEcalRecHitGPU", "EcalRecHitsEE"),
    )

    # SwitchProducer wrapping the legacy ECAL calibrated rechits producer or a converter from SoA to legacy format
    from RecoLocalCalo.EcalRecProducers.ecalRecHitConvertGPU2CPUFormat_cfi import ecalRecHitConvertGPU2CPUFormat as _ecalRecHitConvertGPU2CPUFormat
    process.hltEcalRecHit = SwitchProducerCUDA(
        # legacy producer
        cpu = process.hltEcalRecHit,
        # convert the ECAL calibrated rechits from SoA to legacy format
        cuda = _ecalRecHitConvertGPU2CPUFormat.clone(
            recHitsLabelGPUEB = ("hltEcalRecHitSoA", "EcalRecHitsEB"),
            recHitsLabelGPUEE = ("hltEcalRecHitSoA", "EcalRecHitsEE"),
        )
    )
    """

    # SwitchProducer wrapping the legacy ECAL rechits producer
    # the gpu unpacker does not produce the TPs used for the recovery, so the SwitchProducer alias does not provide them:
    #   - the cpu uncalibrated rechit producer may mark them for recovery, read the TPs explicitly from the legacy unpacker
    #   - the gpu uncalibrated rechit producer does not flag them for recovery, so the TPs are not necessary
    if not isinstance(process.hltEcalRecHit, SwitchProducerCUDA):
        process.hltEcalRecHit = SwitchProducerCUDA(
            cpu = process.hltEcalRecHit.clone(
                triggerPrimitiveDigiCollection = ('hltEcalDigisLegacy', 'EcalTriggerPrimitives')
            )
        )

    elif not hasattr(process.hltEcalRecHit, 'cpu'):
        raise Exception('unsupported configuration: "process.hltEcalRecHit" is a SwitchProducerCUDA, but does not have a "cpu" branch')

    if not hasattr(process.hltEcalRecHit, 'cuda'):
        process.hltEcalRecHit.cuda = process.hltEcalRecHit.cpu.clone(
            triggerPrimitiveDigiCollection = 'unused'
        )


    # Tasks and Sequences

    if 'HLTDoFullUnpackingEgammaEcalWithoutPreshowerTask' not in process.__dict__:
        process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerTask = cms.Task(
            process.hltEcalDigisGPU,             # unpack ECAL digis on gpu
            process.hltEcalDigisLegacy,          # legacy producer, referenced in the SwitchProducer
            process.hltEcalDigis,                # SwitchProducer
            process.hltEcalUncalibRecHitGPU,     # run ECAL local reconstruction and multifit on gpu
            process.hltEcalUncalibRecHitSoA,     # needed by hltEcalPhiSymFilter - copy to host
            process.hltEcalUncalibRecHit,        # needed by hltEcalPhiSymFilter - convert to legacy format
#           process.hltEcalRecHitGPU,            # make ECAL calibrated rechits on gpu
#           process.hltEcalRecHitSoA,            # copy to host
            process.hltEcalDetIdToBeRecovered,   # legacy producer
            process.hltEcalRecHit,               # legacy producer
        )

    elif not isinstance(process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerTask, cms.Task):
        raise Exception('unsupported configuration: "process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerTask" already exists, but it is not a Task')

    # redefine HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence (it was emptied at the start of this function)
    process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence = cms.Sequence(
        process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerTask
    )

    if hasHLTEcalPreshowerSeq:

        if 'HLTPreshowerTask' not in process.__dict__:
            process.HLTPreshowerTask = cms.Task(
                process.hltEcalPreshowerDigis,   # unpack ECAL preshower digis on the host
                process.hltEcalPreshowerRecHit,  # build ECAL preshower rechits on the host
            )

        elif not isinstance(process.HLTPreshowerTask, cms.Task):
            raise Exception('unsupported configuration: "process.HLTPreshowerTask" already exists, but it is not a Task')

        # redefine HLTPreshowerSequence (it was emptied at the start of this function)
        process.HLTPreshowerSequence = cms.Sequence(process.HLTPreshowerTask)

        if 'HLTDoFullUnpackingEgammaEcalTask' not in process.__dict__:
            process.HLTDoFullUnpackingEgammaEcalTask = cms.Task(
                process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerTask,
                process.HLTPreshowerTask,
            )

        elif not isinstance(process.HLTDoFullUnpackingEgammaEcalTask, cms.Task):
            raise Exception('unsupported configuration: "process.HLTDoFullUnpackingEgammaEcalTask" already exists, but it is not a Task')

        # redefine sequences (they were emptied at the start of this function)
        process.HLTDoFullUnpackingEgammaEcalSequence = cms.Sequence(process.HLTDoFullUnpackingEgammaEcalTask)
        process.HLTDoFullUnpackingEgammaEcalMFSequence = cms.Sequence(process.HLTDoFullUnpackingEgammaEcalTask)


    # done
    return process

# customisation for offloading the HCAL local reconstruction via CUDA if a supported gpu is present
def customiseHcalLocalReconstruction(process):

    hasHLTDoLocalHcalSeq = 'HLTDoLocalHcalSequence' in process.__dict__
    if not (hasHLTDoLocalHcalSeq or 'HLTStoppedHSCPLocalHcalReco' in process.__dict__):
        return process

    # FIXME replace the Sequences with empty ones to avoid expanding them during the (re)definition of Modules and EDAliases
    if hasHLTDoLocalHcalSeq:
        process.HLTDoLocalHcalSequence = cms.Sequence()
    process.HLTStoppedHSCPLocalHcalReco = cms.Sequence()


    # Event Setup
    if 'hcalElectronicsMappingGPUESProducer' not in process.__dict__:
        process.load("EventFilter.HcalRawToDigi.hcalElectronicsMappingGPUESProducer_cfi")

    if 'hcalChannelQualityGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalChannelQualityGPUESProducer_cfi")

    if 'hcalGainsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalGainsGPUESProducer_cfi")

    if 'hcalGainWidthsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalGainWidthsGPUESProducer_cfi")

    if 'hcalLUTCorrsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalLUTCorrsGPUESProducer_cfi")

    if 'hcalConvertedPedestalsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalConvertedPedestalsGPUESProducer_cfi")

    if 'hcalConvertedEffectivePedestalsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalConvertedEffectivePedestalsGPUESProducer_cfi")
        process.hcalConvertedEffectivePedestalsGPUESProducer.label0 = "withTopoEff"

    if 'hcalConvertedPedestalWidthsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalConvertedPedestalWidthsGPUESProducer_cfi")

    if 'hcalConvertedEffectivePedestalWidthsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalConvertedEffectivePedestalWidthsGPUESProducer_cfi")
        process.hcalConvertedEffectivePedestalWidthsGPUESProducer.label0 = "withTopoEff"
        process.hcalConvertedEffectivePedestalWidthsGPUESProducer.label1 = "withTopoEff"

    if 'hcalQIECodersGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalQIECodersGPUESProducer_cfi")

    if 'hcalRecoParamsWithPulseShapesGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalRecoParamsWithPulseShapesGPUESProducer_cfi")

    if 'hcalRespCorrsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalRespCorrsGPUESProducer_cfi")

    if 'hcalTimeCorrsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalTimeCorrsGPUESProducer_cfi")

    if 'hcalQIETypesGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalQIETypesGPUESProducer_cfi")

    if 'hcalSiPMParametersGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalSiPMParametersGPUESProducer_cfi")

    if 'hcalSiPMCharacteristicsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalSiPMCharacteristicsGPUESProducer_cfi")

    if 'hcalMahiPulseOffsetsGPUESProducer' not in process.__dict__:
        process.load("RecoLocalCalo.HcalRecProducers.hcalMahiPulseOffsetsGPUESProducer_cfi")


    # Modules and EDAliases

    # The HCAL unpacker running on the gpu supports only the HB and HE digis.
    # So, run the legacy unacker on the cpu, then convert the HB and HE digis
    # to SoA format and copy them to the gpu.

    if 'hltHcalDigisGPU' not in process.__dict__:
        from EventFilter.HcalRawToDigi.hcalDigisProducerGPU_cfi import hcalDigisProducerGPU as _hcalDigisProducerGPU
        process.hltHcalDigisGPU = _hcalDigisProducerGPU.clone(
            hbheDigisLabel = "hltHcalDigis",
            qie11DigiLabel = "hltHcalDigis",
            digisLabelF01HE = "",
            digisLabelF5HB = "",
            digisLabelF3HB = ""
        )

    if 'hltHbherecoGPU' not in process.__dict__:
        # run the HCAL local reconstruction (including Method 0 and MAHI) on gpu
        from RecoLocalCalo.HcalRecProducers.hbheRecHitProducerGPU_cfi import hbheRecHitProducerGPU as _hbheRecHitProducerGPU
        process.hltHbherecoGPU = _hbheRecHitProducerGPU.clone(
            digisLabelF01HE = "hltHcalDigisGPU",
            digisLabelF5HB = "hltHcalDigisGPU",
            digisLabelF3HB = "hltHcalDigisGPU",
            recHitsLabelM0HBHE = ""
        )

    if 'hltHbherecoFromGPU' not in process.__dict__:
        # transfer the HCAL rechits to the cpu, and convert them to the legacy format
        from RecoLocalCalo.HcalRecProducers.hcalCPURecHitsProducer_cfi import hcalCPURecHitsProducer as _hcalCPURecHitsProducer
        process.hltHbherecoFromGPU = _hcalCPURecHitsProducer.clone(
            recHitsM0LabelIn = "hltHbherecoGPU",
            recHitsM0LabelOut = "",
            recHitsLegacyLabelOut = ""
        )


    # SwitchProducer between the legacy producer and the copy from gpu with conversion
    if not isinstance(process.hltHbhereco, SwitchProducerCUDA):
        process.hltHbhereco = SwitchProducerCUDA(
            # legacy producer
            cpu = process.hltHbhereco.clone()
        )

    elif not hasattr(process.hltHbhereco, 'cpu'):
        raise Exception('unsupported configuration: "process.hltHbhereco" is a SwitchProducerCUDA, but does not have a "cpu" branch')

    if not hasattr(process.hltHbhereco, 'cuda'):
        # alias to the rechits converted to legacy format
        process.hltHbhereco.cuda = cms.EDAlias(
            hltHbherecoFromGPU = cms.VPSet(
                cms.PSet(type = cms.string("HBHERecHitsSorted"))
            )
        )


    # Tasks and Sequences

    if hasHLTDoLocalHcalSeq:
        if 'HLTDoLocalHcalTask' not in process.__dict__:
            process.HLTDoLocalHcalTask = cms.Task(
                process.hltHcalDigis,       # legacy producer, unpack HCAL digis on cpu
                process.hltHcalDigisGPU,    # copy to gpu and convert to SoA format
                process.hltHbherecoGPU,     # run the HCAL local reconstruction (including Method 0 and MAHI) on gpu
                process.hltHbherecoFromGPU, # transfer the HCAL rechits to the cpu, and convert them to the legacy format
                process.hltHbhereco,        # SwitchProducer between the legacy producer and the copy from gpu with conversion
                process.hltHfprereco,       # legacy producer
                process.hltHfreco,          # legacy producer
                process.hltHoreco,          # legacy producer
            )

        elif not isinstance(process.HLTDoLocalHcalTask, cms.Task):
            raise Exception('unsupported configuration: "process.HLTDoLocalHcalTask" already exists, but it is not a Task')

        # redefine HLTDoLocalHcalSequence (it was emptied at the start of this function)
        process.HLTDoLocalHcalSequence = cms.Sequence(process.HLTDoLocalHcalTask)

    if 'HLTStoppedHSCPLocalHcalRecoTask' not in process.__dict__:
        process.HLTStoppedHSCPLocalHcalRecoTask = cms.Task(
            process.hltHcalDigis,           # legacy producer, unpack HCAL digis on cpu
            process.hltHcalDigisGPU,        # copy to gpu and convert to SoA format
            process.hltHbherecoGPU,         # run the HCAL local reconstruction (including Method 0 and MAHI) on gpu
            process.hltHbherecoFromGPU,     # transfer the HCAL rechits to the cpu, and convert them to the legacy format
            process.hltHbhereco,            # SwitchProducer between the legacy producer and the copy from gpu with conversion
        )

    elif not isinstance(process.HLTStoppedHSCPLocalHcalRecoTask, cms.Task):
        raise Exception('unsupported configuration: "process.HLTStoppedHSCPLocalHcalRecoTask" already exists, but it is not a Task')

    # redefine HLTStoppedHSCPLocalHcalReco (it was emptied at the start of this function)
    process.HLTStoppedHSCPLocalHcalReco = cms.Sequence(process.HLTStoppedHSCPLocalHcalRecoTask)


    # done
    return process


# customisation to enable pixel triplets instead of quadruplets
def enablePatatrackPixelTriplets(process):

  if 'hltPixelTracksCUDA' in process.__dict__:
      # configure GPU pixel tracks for triplets
      process.hltPixelTracksCUDA.minHitsPerNtuplet = 3
      process.hltPixelTracksCUDA.includeJumpingForwardDoublets = True

  if 'hltPixelTracksSoA' in process.__dict__:
      # configure CPU pixel tracks for triplets
      process.hltPixelTracksSoA.cpu.minHitsPerNtuplet = 3
      process.hltPixelTracksSoA.cpu.includeJumpingForwardDoublets = True

  # done
  return process


# customisation for running the Patatrack reconstruction, with automatic offload via CUDA when a supported gpu is available
def customizeHLTforPatatrack(process):
    process = customiseCommon(process)
    process = customisePixelLocalReconstruction(process)
    process = customisePixelTrackReconstruction(process)
    process = customiseEcalLocalReconstruction(process)
    process = customiseHcalLocalReconstruction(process)
    return process


# customisation for running the Patatrack triplets reconstruction, with automatic offload via CUDA when a supported gpu is available
def customizeHLTforPatatrackTriplets(process):
    process = customiseCommon(process)
    process = customisePixelLocalReconstruction(process)
    process = customisePixelTrackReconstruction(process)
    process = customiseEcalLocalReconstruction(process)
    process = customiseHcalLocalReconstruction(process)
    process = enablePatatrackPixelTriplets(process)
    return process


def _addConsumerPath(process):
    # add to a path all consumers and the tasks that define the producers
    process.Consumer = cms.Path(
        process.HLTBeginSequence +
        process.hltPixelConsumer +
        process.hltEcalConsumer +
        process.hltHbheConsumer,
        process.HLTDoLocalPixelTask,
        process.HLTRecoPixelTracksTask,
        process.HLTRecopixelvertexingTask,
        process.HLTDoFullUnpackingEgammaEcalTask,
        process.HLTDoLocalHcalTask,
    )

    if process.schedule is not None:
        process.schedule.append(process.Consumer)

    # done
    return process


def consumeGPUSoAProducts(process):
    # consume the Pixel tracks and vertices on the GPU in SoA format
    process.hltPixelConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltPixelTracksCUDA', 'hltPixelVerticesCUDA' )
    )

    # consume the ECAL uncalibrated rechits on the GPU in SoA format
    process.hltEcalConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltEcalUncalibRecHitGPU' )
    )

    # consume the HCAL rechits on the GPU in SoA format
    process.hltHbheConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltHbherecoGPU' )
    )

    # add to a path all consumers and the tasks that define the producers
    process = _addConsumerPath(process)

    # done
    return process


def consumeCPUSoAProducts(process):
    # consume the Pixel tracks and vertices on the CPU in SoA format
    process.hltPixelConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltPixelTracksSoA', 'hltPixelVerticesSoA' )
    )

    # consume the ECAL uncalibrated rechits on the CPU in SoA format
    process.hltEcalConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltEcalUncalibRecHitSoA' )
    )

    # consume the HCAL rechits on the CPU in legacy format
    process.hltHbheConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltHbhereco' )
    )

    # add to a path all consumers and the tasks that define the producers
    process = _addConsumerPath(process)

    # done
    return process

def consumeCPULegacyProducts(process):
    # consume the Pixel tracks and vertices on the CPU in legacy format
    process.hltPixelConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltPixelTracks', 'hltPixelVertices' )
    )

    # consume the ECAL runcalibrated echits on the CPU in legacy format
    process.hltEcalConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltEcalUncalibRecHit' )
    )

    # consume the HCAL rechits on the CPU in legacy format
    process.hltHbheConsumer = cms.EDAnalyzer("GenericConsumer",
        eventProducts = cms.untracked.vstring( 'hltHbhereco' )
    )

    # add to a path all consumers and the tasks that define the producers
    process = _addConsumerPath(process)

    # done
    return process
