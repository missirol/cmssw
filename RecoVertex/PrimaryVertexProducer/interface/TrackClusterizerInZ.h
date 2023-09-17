#ifndef TrackClusterizerInZ_h
#define TrackClusterizerInZ_h

/**\class TrackClusterizerInZ 
 
  Description: interface/base class for track clusterizers that separate event tracks into clusters along the beam line

*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

class TrackClusterizerInZ {
public:
  TrackClusterizerInZ() {}
  TrackClusterizerInZ(const edm::ParameterSet& conf) {}
  virtual ~TrackClusterizerInZ(){};

  virtual std::vector<std::vector<reco::TransientTrack> > clusterize(const std::vector<reco::TransientTrack>& tracks) const = 0;
  virtual std::vector<TransientVertex> vertices(const std::vector<reco::TransientTrack>& tracks) const = 0;
};

#endif
