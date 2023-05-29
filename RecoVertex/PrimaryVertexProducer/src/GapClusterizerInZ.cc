#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"

namespace {

  bool recTrackLessZ(const reco::TransientTrack& tk1, const reco::TransientTrack& tk2) {
    return tk1.stateAtBeamLine().trackStateAtPCA().position().z() <
           tk2.stateAtBeamLine().trackStateAtPCA().position().z();
  }

}  // namespace

GapClusterizerInZ::GapClusterizerInZ(const edm::ParameterSet& conf):
  zSep_{(float) conf.getParameter<double>("zSeparation")},
  verbose_{conf.getUntrackedParameter<bool>("verbose")} {
  if (verbose_) {
    std::cout << "TrackClusterizerInZ:  algorithm=gap, zSeparation=" << zSep_ << std::endl;
  }
}

std::vector<std::vector<reco::TransientTrack>> GapClusterizerInZ::clusterize(const std::vector<reco::TransientTrack>& tracks) const {
  std::vector<reco::TransientTrack> tks = tracks;  // copy to be sorted

  std::vector<std::vector<reco::TransientTrack>> clusters;
  if (tks.empty())
    return clusters;

  // sort in increasing order of z
  stable_sort(tks.begin(), tks.end(), recTrackLessZ);

  // init first cluster
  std::vector<reco::TransientTrack>::const_iterator it = tks.begin();
  std::vector<reco::TransientTrack> currentCluster;
  currentCluster.push_back(*it);

  it++;
  for (; it != tks.end(); it++) {
    double zPrev = currentCluster.back().stateAtBeamLine().trackStateAtPCA().position().z();
    double zCurr = (*it).stateAtBeamLine().trackStateAtPCA().position().z();

    if (std::abs(zCurr - zPrev) < zSeparation()) {
      // close enough ? cluster together
      currentCluster.push_back(*it);
    } else {
      // store current cluster, start new one
      clusters.push_back(currentCluster);
      currentCluster.clear();
      currentCluster.push_back(*it);
    }
  }

  // store last cluster
  clusters.push_back(currentCluster);

  return clusters;
}

std::vector<TransientVertex> GapClusterizerInZ::vertices(const std::vector<reco::TransientTrack>& tracks) const {
  /* repackage track clusters, compatibility with newer clusterizers */
  std::vector<TransientVertex> primary_vertices;
  auto trackClusters = clusterize(tracks);

  GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);
  for (auto& vertexTracks : trackClusters) {
    GlobalPoint position(0, 0, 0);  // dummy
    primary_vertices.push_back(TransientVertex(position, dummyError, vertexTracks, 0));
  }

  return primary_vertices;
}

void GapClusterizerInZ::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<double>("zSeparation", 1.0);
  desc.addUntracked<bool>("verbose", false);
}
