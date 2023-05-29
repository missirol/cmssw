#ifndef RecoVertex_PrimaryVertexProducer_VertexTimeAlgorithmBase_h
#define RecoVertex_PrimaryVertexProducer_VertexTimeAlgorithmBase_h

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
  class ParameterSetDescription;
  class ConsumesCollector;
}  // namespace edm

class VertexTimeAlgorithmBase {
public:
  VertexTimeAlgorithmBase(edm::ParameterSet const& iConfig, edm::ConsumesCollector& iCC) {}
  VertexTimeAlgorithmBase(edm::ParameterSet const& iConfig, edm::ConsumesCollector&& iCC) : VertexTimeAlgorithmBase(iConfig, iCC) {}
  VertexTimeAlgorithmBase(VertexTimeAlgorithmBase const&) = delete;
  VertexTimeAlgorithmBase& operator=(VertexTimeAlgorithmBase const&) = delete;
  virtual ~VertexTimeAlgorithmBase() = default;

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {}

  virtual void setEvent(edm::Event& iEvent, edm::EventSetup const& iSetup) = 0;

  /**
   * estimate the vertex time and time uncertainty for transient vertex
   * 
   * returns true when a valid time has been determined, otherwise return false
   */
  virtual bool vertexTime(float& vtxTime, float& vtxTimeError, TransientVertex const& vtx) const = 0;

  /**
   * replace the vertices in the input vector by new vertices with time coordinates
   * determined by the vertexTime method
   * this implementation does not alter the weights from the previous fit
   * must be overridden to change weights, coordinates, tracklists or to add or remove vertices
   */
  virtual void fill_vertex_times(std::vector<TransientVertex>& pvs) const;
};

#endif
