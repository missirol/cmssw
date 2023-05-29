#ifndef RecoVertex_PrimaryVertexProducer_TrackFilterForPVFinding_h
#define RecoVertex_PrimaryVertexProducer_TrackFilterForPVFinding_h

/**\class TrackFilterForPVFinding 
 
  Description: track selection for PV finding

*/
#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFindingBase.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}

class TrackFilterForPVFinding : public TrackFilterForPVFindingBase {
public:
  TrackFilterForPVFinding(const edm::ParameterSet& conf);

  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  bool operator()(const reco::TransientTrack& tracks) const;
  std::vector<reco::TransientTrack> select(const std::vector<reco::TransientTrack>& tracks) const override;
  std::vector<reco::TransientTrack> selectTight(const std::vector<reco::TransientTrack>& tracks,
                                                double minPtTight) const;

private:
  float maxD0Sig_, minPt_, maxEta_;
  float maxD0Error_, maxDzError_;
  int minSiLayers_, minPxLayers_;
  float maxNormChi2_;
  reco::TrackBase::TrackQuality quality_;
};

#endif
