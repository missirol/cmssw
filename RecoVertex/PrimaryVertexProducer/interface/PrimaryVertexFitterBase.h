#ifndef PrimaryVertexFitterBase_h
#define PrimaryVertexFitterBase_h

/**\class PrimaryVertexFitterBase
 
  Description: base class for primary vertex fitters

*/

//#include "FWCore/ParameterSet/interface/ParameterSet.h"
namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}

namespace reco {
  class BeamSpot;
  class TransientTrack;
}
class TransientVertex;

//#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class PrimaryVertexFitterBase {
 public:
  PrimaryVertexFitterBase(const edm::ParameterSet &conf){}
  PrimaryVertexFitterBase(){}
  virtual ~PrimaryVertexFitterBase() = default;
  //virtual std::vector<TransientVertex> vertices(const std::vector<TransientVertex> &, const reco::BeamSpot &, const bool ) = 0;
  virtual std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack> &, const std::vector<TransientVertex> &, const reco::BeamSpot &, const bool ) = 0;
};
#endif
