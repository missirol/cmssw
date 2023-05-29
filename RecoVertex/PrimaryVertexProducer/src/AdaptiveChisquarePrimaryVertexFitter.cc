#include "RecoVertex/PrimaryVertexProducer/interface/AdaptiveChisquarePrimaryVertexFitter.h"

std::vector<TransientVertex> AdaptiveChisquarePrimaryVertexFitter::fit(std::vector<reco::TransientTrack> const&, std::vector<TransientVertex> const& clusters, reco::BeamSpot const& beamspot, bool const useBeamConstraint) {
  // fit the clusters one-by-one
  std::vector<TransientVertex> pvs;
  pvs.reserve(clusters.size());

  for (auto const& cluster : clusters) {
    auto const& transTracks = cluster.originalTracks();

    if (transTracks.size() <= 1) {
      continue;
    }

    auto const result = fitter_.fit(transTracks, {cluster}, beamspot, useBeamConstraint);

    if (result.empty() or not result[0].isValid()) {
      continue;
    }

    pvs.emplace_back(result[0]);
  }

  return pvs;
}
