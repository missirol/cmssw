#ifndef RecoVertex_PrimaryVertexProducer_SequentialPrimaryVertexFitterAdapter_h
#define RecoVertex_PrimaryVertexProducer_SequentialPrimaryVertexFitterAdapter_h

/**\class SequentialPrimaryVertexFitterAdapter
 
  Description: Adapter class for Kalman and Adaptive vertex fitters 

*/

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"

class SequentialPrimaryVertexFitterAdapter : public PrimaryVertexFitterBase {
public:
  SequentialPrimaryVertexFitterAdapter() : fitter(nullptr) {}
  SequentialPrimaryVertexFitterAdapter(VertexFitter<5>* vertex_fitter) : fitter(vertex_fitter) {}
  ~SequentialPrimaryVertexFitterAdapter() override = default;

  std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack>&,
                                   const std::vector<TransientVertex>& clusters,
                                   const reco::BeamSpot& beamspot,
                                   const bool useBeamConstraint) override {
    std::vector<TransientVertex> pvs;
    for (auto const& cluster : clusters) {
      auto const& tracklist = cluster.originalTracks();

      if (tracklist.size() <= 1) {
        continue;
      }

      auto const& v = useBeamConstraint ? fitter->vertex(tracklist, beamspot) : fitter->vertex(tracklist);

      if (v.isValid()) {
        pvs.emplace_back(v);
      }
    }
    return pvs;
  };

protected:
  // configuration
  VertexFitter<5>* fitter;  // Kalman or Adaptive
};

#endif
