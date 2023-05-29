#ifndef RecoVertex_PrimaryVertexProducer_AdaptiveChisquarePrimaryVertexFitter_h
#define RecoVertex_PrimaryVertexProducer_AdaptiveChisquarePrimaryVertexFitter_h

/**\class AdaptiveChisquarePrimaryVertexFitter
 
  Description: Adapter for using the MultiPrimaryVertexFitter to fit vertices one-by-one instead of all simultaneously

*/
#include <memory>
#include <vector>

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/MultiPrimaryVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

class AdaptiveChisquarePrimaryVertexFitter : public PrimaryVertexFitterBase {
public:
  AdaptiveChisquarePrimaryVertexFitter(double const chi2cutoff, double const mintrkweight)
    : fitter_(chi2cutoff, mintrkweight) {}

  std::vector<TransientVertex> fit(std::vector<reco::TransientTrack> const&,
                                   std::vector<TransientVertex> const& clusters,
                                   reco::BeamSpot const& beamspot,
                                   bool const useBeamConstraint) override;
protected:
  MultiPrimaryVertexFitter fitter_;
};
#endif
