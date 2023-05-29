#ifndef RecoVertex_PrimaryVertexProducer_SequentialPrimaryVertexFitterAdapter_h
#define RecoVertex_PrimaryVertexProducer_SequentialPrimaryVertexFitterAdapter_h

/**\class SequentialPrimaryVertexFitterAdapter
 
  Description: Adapter class for Kalman and Adaptive vertex fitters 

*/
#include <memory>
#include <vector>

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"

class SequentialPrimaryVertexFitterAdapter : public PrimaryVertexFitterBase {
public:
  SequentialPrimaryVertexFitterAdapter(std::unique_ptr<VertexFitter<5>> fitter) : fitter_(std::move(fitter)) {}
  std::vector<TransientVertex> fit(std::vector<reco::TransientTrack> const&, std::vector<TransientVertex> const& clusters, reco::BeamSpot const& beamspot, bool const useBeamConstraint) override;

protected:
  // vertex fitter (Kalman or Adaptive)
  std::unique_ptr<VertexFitter<5>> fitter_;
};

#endif
