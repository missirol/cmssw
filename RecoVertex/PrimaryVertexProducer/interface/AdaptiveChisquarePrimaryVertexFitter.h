#ifndef AdaptiveChisquarePrimaryVertexFitter_h
#define AdaptiveChisquarePrimaryVertexFitter_h

/**\class AdaptiveChisquarePrimaryVertexFitter
 
  Description: Adapter for using the MultiPrimaryVertexFitter to fit vertices one-by-one instead of all simultaneously

*/

//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/MultiPrimaryVertexFitter.h"

class AdaptiveChisquarePrimaryVertexFitter : public PrimaryVertexFitterBase {
public:
  AdaptiveChisquarePrimaryVertexFitter(const double chi2cutoff, const double mintrkweight)
      : fitter(new MultiPrimaryVertexFitter(chi2cutoff, mintrkweight)){};
  ~AdaptiveChisquarePrimaryVertexFitter() override = default;  // FIXME, delete the fitter

  std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack>& dummy,
                                   const std::vector<TransientVertex>& clusters,
                                   const reco::BeamSpot& beamspot,
                                   const bool useBeamConstraint) override {
    // fit the clusters one-by-one
    std::vector<TransientVertex> pvs;
    std::vector<TransientVertex> seed(1);
    for (auto& cluster : clusters) {
      if ((useBeamConstraint && (cluster.originalTracks().size() > 1)) || (cluster.originalTracks().size() > 1)) {
        const std::vector<reco::TransientTrack>& tracklist = cluster.originalTracks();
        seed[0] = cluster;
        auto result = fitter->fit(tracklist, seed, beamspot, useBeamConstraint);
        if ((!result.empty()) && (result[0].isValid())) {
          pvs.push_back(result[0]);
        }
      }
    }
    return pvs;
  };

protected:
  // configuration
  MultiPrimaryVertexFitter* fitter;  // Kalman or Adaptive
};
#endif
