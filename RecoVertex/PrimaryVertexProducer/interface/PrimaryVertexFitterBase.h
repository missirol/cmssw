#ifndef RecoVertex_PrimaryVertexProducer_PrimaryVertexFitterBase_h
#define RecoVertex_PrimaryVertexProducer_PrimaryVertexFitterBase_h

/**\class PrimaryVertexFitterBase

  Description: base class for primary vertex fitters

*/
namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

namespace reco {
  class BeamSpot;
  class TransientTrack;
}  // namespace reco

class TransientVertex;

class PrimaryVertexFitterBase {
public:
  PrimaryVertexFitterBase() {}
  PrimaryVertexFitterBase(const edm::ParameterSet &conf) {}
  virtual ~PrimaryVertexFitterBase() = default;
  virtual std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack> &,
                                           const std::vector<TransientVertex> &,
                                           const reco::BeamSpot &,
                                           const bool) = 0;
};

#endif
