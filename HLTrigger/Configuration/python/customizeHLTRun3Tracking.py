import copy
import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.common import *
from HLTrigger.Configuration.customizeHLTforPatatrack import *

def customizeHLTRun3Tracking(process):

    for producer in producers_by_type(process, "TrackWithVertexSelector"):
        if not hasattr(producer,'numberOfValidHitsForGood'):
            producer.numberOfValidHitsForGood = cms.uint32(999)
        if not hasattr(producer,'numberOfValidPixelHitsForGood'):
            producer.numberOfValidPixelHitsForGood = cms.uint32(999)
        if not hasattr(producer,'zetaVtxScale'):
            producer.zetaVtxScale = cms.double(1.0)
        if not hasattr(producer,'rhoVtxScale'):
            producer.rhoVtxScale = cms.double(1.0)
        if not hasattr(producer,'zetaVtxSig'):
            producer.zetaVtxSig = cms.double(999.0)
        if not hasattr(producer,'rhoVtxSig'):
            producer.rhoVtxSig = cms.double(999.0)
    
    process = customizeHLTforPatatrack(process)    
    if 'hltPixelTracksCUDA' not in process.__dict__:
        return process
    process.hltPixelTracksCUDA.includeJumpingForwardDoublets = cms.bool(True)
    process.hltPixelTracksCUDA.minHitsPerNtuplet             = cms.uint32(3)
    process.hltPixelTracksCUDA.idealConditions = cms.bool(True)

    if 'hltPixelTracksSoA' not in process.__dict__:
        return process
    process.hltPixelTracksSoA.cpu.includeJumpingForwardDoublets = cms.bool(True)
    process.hltPixelTracksSoA.cpu.minHitsPerNtuplet             = cms.uint32(3)
    process.hltPixelTracksSoA.cpu.idealConditions = cms.bool(True)

    if 'HLTIter0PSetTrajectoryFilterIT' not in process.__dict__:
        return process
    process.HLTIter0PSetTrajectoryFilterIT.minHitsMinPt        = cms.int32(3)
    process.HLTIter0PSetTrajectoryFilterIT.minimumNumberOfHits = cms.int32(3)

    if 'hltSiStripRawToClustersFacility' not in process.__dict__:
        return process
    process.hltSiStripRawToClustersFacility.onDemand = cms.bool( False )
    
    if 'hltIter0PFLowPixelSeedsFromPixelTracks' not in process.__dict__:
        return process
    process.hltIter0PFLowPixelSeedsFromPixelTracks.includeFourthHit = cms.bool(True)

    if 'hltIter0PFlowTrackCutClassifier' not in process.__dict__:
        return process
    del process.hltIter0PFlowTrackCutClassifier
    process.hltIter0PFlowTrackCutClassifier = cms.EDProducer(
        "TrackCutClassifier",
        src = cms.InputTag("hltIter0PFlowCtfWithMaterialTracks"),
        beamspot = cms.InputTag("hltOnlineBeamSpot"),
        vertices = cms.InputTag("hltTrimmedPixelVertices"),
        qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
        mva = cms.PSet(
            minPixelHits = cms.vint32(0, 0, 0),
            maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
            dr_par = cms.PSet(
                d0err = cms.vdouble(0.003, 0.003, 0.003),
                dr_par2 = cms.vdouble(3.40282346639e+38, 0.6, 0.6),
                dr_par1 = cms.vdouble(3.40282346639e+38, 0.8, 0.8),
                dr_exp = cms.vint32(4, 4, 4),
                d0err_par = cms.vdouble(0.001, 0.001, 0.001)
            ),
            maxLostLayers = cms.vint32(1, 1, 1),
            min3DLayers = cms.vint32(0, 0, 0),
            dz_par = cms.PSet(
                dz_par1 = cms.vdouble(3.40282346639e+38, 0.75, 0.75),
                dz_par2 = cms.vdouble(3.40282346639e+38, 0.5, 0.5),
                dz_exp = cms.vint32(4, 4, 4)
            ),
            minNVtxTrk = cms.int32(3),
            maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
            minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
            maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
            maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
            maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
            minLayers = cms.vint32(3, 3, 3)
        ),
        ignoreVertices = cms.bool(False)
    )
        
    if 'hltMergedTracks' not in process.__dict__:
        return process
    del process.hltMergedTracks
    process.hltMergedTracks = process.hltIter0PFlowTrackSelectionHighPurity.clone()

    if 'HLTIterativeTrackingIteration0' not in process.__dict__:
        return process
    process.HLTIterativeTrackingIteration0 = cms.Sequence( 
        process.hltIter0PFLowPixelSeedsFromPixelTracks+
        process.hltIter0PFlowCkfTrackCandidates+
        process.hltIter0PFlowCtfWithMaterialTracks+
        process.hltIter0PFlowTrackCutClassifier+
        process.hltMergedTracks
    )
    
    if 'HLTIterativeTrackingIter02' not in process.__dict__:
        return process
    process.HLTIterativeTrackingIter02 = cms.Sequence( 
        process.HLTIterativeTrackingIteration0 
    )
    
    if 'MC_ReducedIterativeTracking_v12' not in process.__dict__:
        return process
    process.MC_ReducedIterativeTracking_v12 = cms.Path( 
        process.HLTBeginSequence +
        process.hltPreMCReducedIterativeTracking +
        process.HLTDoLocalPixelSequence +
        process.HLTRecopixelvertexingSequence +
        process.HLTDoLocalStripSequence +
        process.HLTIterativeTrackingIter02 +
        process.HLTEndSequence
    )

    return process



def customizeHLTRun3TrackingAllPixelVertices(process):

    for producer in producers_by_type(process, "TrackWithVertexSelector"):
        if not hasattr(producer,'numberOfValidHitsForGood'):
            producer.numberOfValidHitsForGood = cms.uint32(999)
        if not hasattr(producer,'numberOfValidPixelHitsForGood'):
            producer.numberOfValidPixelHitsForGood = cms.uint32(999)
        if not hasattr(producer,'zetaVtxScale'):
            producer.zetaVtxScale = cms.double(1.0)
        if not hasattr(producer,'rhoVtxScale'):
            producer.rhoVtxScale = cms.double(1.0)
        if not hasattr(producer,'zetaVtxSig'):
            producer.zetaVtxSig = cms.double(999.0)
        if not hasattr(producer,'rhoVtxSig'):
            producer.rhoVtxSig = cms.double(999.0)
    
    process = customizeHLTforPatatrack(process)    
    if 'hltPixelTracksCUDA' not in process.__dict__:
        return process
    process.hltPixelTracksCUDA.includeJumpingForwardDoublets = cms.bool(True)
    process.hltPixelTracksCUDA.minHitsPerNtuplet             = cms.uint32(3)
    process.hltPixelTracksCUDA.idealConditions = cms.bool(True)

    if 'hltPixelTracksSoA' not in process.__dict__:
        return process
    process.hltPixelTracksSoA.cpu.includeJumpingForwardDoublets = cms.bool(True)
    process.hltPixelTracksSoA.cpu.minHitsPerNtuplet             = cms.uint32(3)
    process.hltPixelTracksSoA.cpu.idealConditions = cms.bool(True)

    if 'HLTIter0PSetTrajectoryFilterIT' not in process.__dict__:
        return process
    process.HLTIter0PSetTrajectoryFilterIT.minHitsMinPt        = cms.int32(3)
    process.HLTIter0PSetTrajectoryFilterIT.minimumNumberOfHits = cms.int32(3)

    if 'hltSiStripRawToClustersFacility' not in process.__dict__:
        return process
    process.hltSiStripRawToClustersFacility.onDemand = cms.bool( False )
    
    process.hltPixelTracksClean = cms.EDProducer("TrackWithVertexSelector",
        src = cms.InputTag('hltPixelTracks'),
        etaMin = cms.double(0.0),
        etaMax = cms.double(5.0),
        ptMin = cms.double(0.0),
        ptMax = cms.double(999999.),
        d0Max = cms.double(999.),
        dzMax = cms.double(999.),
        normalizedChi2 = cms.double(999999.),
        numberOfValidHits = cms.uint32(0),
        numberOfValidPixelHits = cms.uint32(0),
        numberOfValidHitsForGood = cms.uint32(4),
        numberOfValidPixelHitsForGood = cms.uint32(4),
        numberOfLostHits = cms.uint32(999),
        ptErrorCut = cms.double(999999.),
        quality = cms.string("any"),
        useVtx = cms.bool(True),
        vertexTag = cms.InputTag('hltPixelVertices'),
        timesTag = cms.InputTag(''),
        timeResosTag = cms.InputTag(''),
        nVertices = cms.uint32(100),
        vtxFallback = cms.bool(True),
        zetaVtx = cms.double(0.3),
        rhoVtx = cms.double(0.1),
        zetaVtxScale = cms.double(5.0),
        rhoVtxScale = cms.double(2.5), 
        zetaVtxSig = cms.double(2.0),
        rhoVtxSig = cms.double(2.0),
        nSigmaDtVertex = cms.double(0),
        copyExtras = cms.untracked.bool(True),
        copyTrajectories = cms.untracked.bool(False)
    )

    if 'hltIter0PFLowPixelSeedsFromPixelTracks' not in process.__dict__:
        return process
    process.hltIter0PFLowPixelSeedsFromPixelTracks.InputVertexCollection = cms.InputTag( "hltPixelVertices" )
    process.hltIter0PFLowPixelSeedsFromPixelTracks.InputCollection = cms.InputTag("hltPixelTracksClean")
    process.hltIter0PFLowPixelSeedsFromPixelTracks.includeFourthHit = cms.bool(True)

    if 'hltIter0PFlowTrackCutClassifier' not in process.__dict__:
        return process
    del process.hltIter0PFlowTrackCutClassifier
    process.hltIter0PFlowTrackCutClassifier = cms.EDProducer(
        "TrackCutClassifier",
        src = cms.InputTag("hltIter0PFlowCtfWithMaterialTracks"),
        beamspot = cms.InputTag("hltOnlineBeamSpot"),
        vertices = cms.InputTag("hltTrimmedPixelVertices"),
        qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
        mva = cms.PSet(
            minPixelHits = cms.vint32(0, 0, 0),
            maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
            dr_par = cms.PSet(
                d0err = cms.vdouble(0.003, 0.003, 0.003),
                dr_par2 = cms.vdouble(3.40282346639e+38, 0.6, 0.6),
                dr_par1 = cms.vdouble(3.40282346639e+38, 0.8, 0.8),
                dr_exp = cms.vint32(4, 4, 4),
                d0err_par = cms.vdouble(0.001, 0.001, 0.001)
            ),
            maxLostLayers = cms.vint32(1, 1, 1),
            min3DLayers = cms.vint32(0, 0, 0),
            dz_par = cms.PSet(
                dz_par1 = cms.vdouble(3.40282346639e+38, 0.75, 0.75),
                dz_par2 = cms.vdouble(3.40282346639e+38, 0.5, 0.5),
                dz_exp = cms.vint32(4, 4, 4)
            ),
            minNVtxTrk = cms.int32(3),
            maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
            minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
            maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
            maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
            maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
            minLayers = cms.vint32(3, 3, 3)
        ),
        ignoreVertices = cms.bool(False)
    )
    process.hltIter0PFlowTrackCutClassifier.vertices = cms.InputTag("hltPixelVertices")
        
    if 'hltMergedTracks' not in process.__dict__:
        return process
    del process.hltMergedTracks
    process.hltMergedTracks = process.hltIter0PFlowTrackSelectionHighPurity.clone()

    if 'HLTIterativeTrackingIteration0' not in process.__dict__:
        return process
    process.HLTIterativeTrackingIteration0 = cms.Sequence( 
        process.hltPixelTracksClean+
        process.hltIter0PFLowPixelSeedsFromPixelTracks+
        process.hltIter0PFlowCkfTrackCandidates+
        process.hltIter0PFlowCtfWithMaterialTracks+
        process.hltIter0PFlowTrackCutClassifier+
        process.hltMergedTracks
    )
    
    if 'HLTIterativeTrackingIter02' not in process.__dict__:
        return process
    process.HLTIterativeTrackingIter02 = cms.Sequence( 
        process.HLTIterativeTrackingIteration0 
    )
    
    if 'MC_ReducedIterativeTracking_v12' not in process.__dict__:
        return process
    process.MC_ReducedIterativeTracking_v12 = cms.Path( 
        process.HLTBeginSequence +
        process.hltPreMCReducedIterativeTracking +
        process.HLTDoLocalPixelSequence +
        process.HLTRecopixelvertexingSequence +
        process.HLTDoLocalStripSequence +
        process.HLTIterativeTrackingIter02 +
        process.HLTEndSequence
    )

    return process
