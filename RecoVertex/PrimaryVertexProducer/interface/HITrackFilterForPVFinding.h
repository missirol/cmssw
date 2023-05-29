#ifndef RecoVertex_PrimaryVertexProducer_HITrackFilterForPVFinding_h
#define RecoVertex_PrimaryVertexProducer_HITrackFilterForPVFinding_h

/**\class HITrackFilterForPVFinding 
 
  Description: selects tracks for primary vertex reconstruction using th TrackFilterForPVFinding,
  returns the input set of tracks if less than NumTracksThreshold tracks were selected

*/
#include <limits>

#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"

namespace edm {
  class ParameterSet;
}

class HITrackFilterForPVFinding : public TrackFilterForPVFinding {
public:
  HITrackFilterForPVFinding(const edm::ParameterSet& conf) : TrackFilterForPVFinding(conf),
    NumTracksThreshold_{(unsigned int) conf.getParameter<int>("numTracksThreshold")},
    MaxNumTracksThreshold_{(unsigned int) conf.getParameter<int>("maxNumTracksThreshold")},
    minPtTight_{conf.getParameter<double>("minPtTight")} {}

  // override the select method
  std::vector<reco::TransientTrack> select(const std::vector<reco::TransientTrack>& tracks) const override {
    std::vector<reco::TransientTrack> seltks = TrackFilterForPVFinding::select(tracks);
    if (seltks.size() < NumTracksThreshold_) {
      return tracks;
    } else if (seltks.size() > MaxNumTracksThreshold_) {
      std::vector<reco::TransientTrack> seltksTight = TrackFilterForPVFinding::selectTight(tracks, minPtTight_);
      if (seltksTight.size() >= NumTracksThreshold_)
        return seltksTight;
    }

    return seltks;
  }

  static void fillPSetDescription(edm::ParameterSetDescription& desc) {
    TrackFilterForPVFinding::fillPSetDescription(desc);
    desc.add<int>("numTracksThreshold", 0);  // HI only
    desc.add<int>("maxNumTracksThreshold", std::numeric_limits<int>::max());
    desc.add<double>("minPtTight", 0.0);
  }

private:
  unsigned int const NumTracksThreshold_;
  unsigned int const MaxNumTracksThreshold_;
  double const minPtTight_;
};

#endif
