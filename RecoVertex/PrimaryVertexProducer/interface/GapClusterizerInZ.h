#ifndef RecoVertex_PrimaryVertexProducer_GapClusterizerInZ_h
#define RecoVertex_PrimaryVertexProducer_GapClusterizerInZ_h

/**\class GapClusterizerInZ
 
  Description: separates event tracks into clusters along the beam line

*/
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

class GapClusterizerInZ : public TrackClusterizerInZ {
public:
  GapClusterizerInZ(const edm::ParameterSet& conf);
  ~GapClusterizerInZ() override = default;

  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  std::vector<std::vector<reco::TransientTrack> > clusterize(
      const std::vector<reco::TransientTrack>& tracks) const override;

  float zSeparation() const { return zSep_; }

  std::vector<TransientVertex> vertices(const std::vector<reco::TransientTrack>& tracks) const override;

private:
  float const zSep_;
  bool const verbose_;
};

#endif
