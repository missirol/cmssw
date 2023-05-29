#ifndef RecoVertex_PrimaryVertexProducer_TrackClusterizerInZ_h
#define RecoVertex_PrimaryVertexProducer_TrackClusterizerInZ_h

/**\class TrackClusterizerInZ 
 
  Description: interface/base class for track clusterizers that separate event tracks into clusters along the beam line

*/
#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class TrackClusterizerInZ {
public:
  TrackClusterizerInZ() {}
  TrackClusterizerInZ(const edm::ParameterSet& conf) {}
  virtual ~TrackClusterizerInZ() = default;

  virtual std::vector<TransientVertex> vertices(const std::vector<reco::TransientTrack>& tracks) const = 0;
  virtual std::vector<std::vector<reco::TransientTrack> > clusterize(
      const std::vector<reco::TransientTrack>& tracks) const = 0;
};

#endif
