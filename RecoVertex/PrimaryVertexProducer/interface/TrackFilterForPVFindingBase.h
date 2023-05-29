#ifndef RecoVertex_PrimaryVertexProducer_TrackFilterForPVFindingBase_h
#define RecoVertex_PrimaryVertexProducer_TrackFilterForPVFindingBase_h

/**\class TrackFilterForPVFindingBase
 
  Description: base class for track selection

*/
#include <vector>

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TrackFilterForPVFindingBase {
public:
  TrackFilterForPVFindingBase() {}
  TrackFilterForPVFindingBase(const edm::ParameterSet& conf) {}
  virtual ~TrackFilterForPVFindingBase() = default;

  virtual std::vector<reco::TransientTrack> select(const std::vector<reco::TransientTrack>& tracks) const = 0;
};

#endif
