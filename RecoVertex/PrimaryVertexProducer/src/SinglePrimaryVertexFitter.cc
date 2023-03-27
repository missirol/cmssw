#include "RecoVertex/PrimaryVertexProducer/interface/SinglePrimaryVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

//SinglePrimaryVertexFitter::SinglePrimaryVertexFitter(){};


std::vector<TransientVertex> SinglePrimaryVertexFitter::fit(const std::vector<reco::TransientTrack> & tracks, const std::vector<TransientVertex> & clusters, const reco::BeamSpot & beamspot, const bool useBeamConstraint){

  std::cout << "SinglePrimaryVertexFitter::fit  clusters.size=" << clusters.size() << "  use beamspot = " << useBeamConstraint << std::endl;
  
  std::vector<TransientVertex> pvs;
  for(auto & cluster : clusters){
    std::vector<reco::TransientTrack> tracklist = cluster.originalTracks();
    TransientVertex v;
    if (useBeamConstraint && (tracklist.size() > 1)) {
      v = fitter->vertex(tracklist, beamspot);
    } else if (!(useBeamConstraint) && (tracklist.size() > 1)) {
      v = fitter->vertex(tracklist);
    }  // else: no fit ==> v.isValid()=False
    
    if (v.isValid()){
      pvs.push_back(v);
    }
  }
  return pvs;
}
