#ifndef SinglePrimaryVertexFitter_h
#define SinglePrimaryVertexFitter_h

/**\class SinglePrimaryVertexFitter
 
  Description: Adapter class for Kalman and Adaptive vertex fitters

*/

//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"

class SinglePrimaryVertexFitter : public PrimaryVertexFitterBase {
 public:
 SinglePrimaryVertexFitter():fitter(nullptr){};
 SinglePrimaryVertexFitter(VertexFitter<5>* vertex_fitter):fitter(vertex_fitter){};
   ~SinglePrimaryVertexFitter() override = default;

  //std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack> &, const std::vector<TransientVertex> &, const reco::BeamSpot &, const bool );
  std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack> &, const std::vector<TransientVertex> &, const reco::BeamSpot &, const bool) override;

 protected:
  // configuration
  VertexFitter<5>* fitter;  // Kalman or Adaptive
};
#endif
