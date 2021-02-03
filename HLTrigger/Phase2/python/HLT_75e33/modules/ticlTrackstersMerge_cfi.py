import FWCore.ParameterSet.Config as cms

ticlTrackstersMerge = cms.EDProducer("TrackstersMergeProducer",
    cosangle_align = cms.double(0.9945),
    debug = cms.bool(True),
    e_over_h_threshold = cms.double(1),
    eid_graph_path = cms.string('RecoHGCal/TICL/data/tf_models/energy_id_v0.pb'),
    eid_input_name = cms.string('input'),
    eid_min_cluster_energy = cms.double(1),
    eid_n_clusters = cms.int32(10),
    eid_n_layers = cms.int32(50),
    eid_output_name_energy = cms.string('output/regressed_energy'),
    eid_output_name_id = cms.string('output/id_probabilities'),
    eta_bin_window = cms.int32(1),
    halo_max_distance2 = cms.double(4),
    layer_clusters = cms.InputTag("hgcalLayerClusters"),
    mightGet = cms.optional.untracked.vstring,
    optimiseAcrossTracksters = cms.bool(True),
    phi_bin_window = cms.int32(1),
    pt_neutral_threshold = cms.double(2),
    pt_sigma_high = cms.double(2),
    pt_sigma_low = cms.double(2),
    resol_calo_offset_em = cms.double(1.5),
    resol_calo_offset_had = cms.double(1.5),
    resol_calo_scale_em = cms.double(0.15),
    resol_calo_scale_had = cms.double(0.15),
    seedingTrk = cms.InputTag("ticlSeedingTrk"),
    track_max_eta = cms.double(3),
    track_max_missing_outerhits = cms.int32(5),
    track_min_eta = cms.double(1.48),
    track_min_pt = cms.double(1),
    tracks = cms.InputTag("generalTracks"),
    trackstersem = cms.InputTag("ticlTrackstersEM"),
    trackstershad = cms.InputTag("ticlTrackstersHAD"),
    tracksterstrk = cms.InputTag("ticlTrackstersTrk"),
    tracksterstrkem = cms.InputTag("ticlTrackstersTrkEM")
)
