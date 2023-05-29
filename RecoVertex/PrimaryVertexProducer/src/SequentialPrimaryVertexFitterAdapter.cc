#include "RecoVertex/PrimaryVertexProducer/interface/SequentialPrimaryVertexFitterAdapter.h"

std::vector<TransientVertex> SequentialPrimaryVertexFitterAdapter::fit(std::vector<reco::TransientTrack> const&, std::vector<TransientVertex> const& clusters, reco::BeamSpot const& beamspot, bool const useBeamConstraint) {
  std::vector<TransientVertex> pvs;
  pvs.reserve(clusters.size());

  for (auto const& cluster : clusters) {
    auto const& transTracks = cluster.originalTracks();

    if (transTracks.size() <= 1) {
      continue;
    }

    auto const&& vtx = useBeamConstraint ? fitter_->vertex(transTracks, beamspot) : fitter_->vertex(transTracks);

    if (vtx.isValid()) {
      pvs.emplace_back(vtx);
    }
  }

  return pvs;
};
