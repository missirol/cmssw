#ifndef RecoVertex_PrimaryVertexProducer_VertexTimeAlgorithmLegacy4D_h
#define RecoVertex_PrimaryVertexProducer_VertexTimeAlgorithmLegacy4D_h

#include "VertexTimeAlgorithmBase.h"

class VertexTimeAlgorithmLegacy4D : public VertexTimeAlgorithmBase {
public:
  VertexTimeAlgorithmLegacy4D(edm::ParameterSet const& iConfig, edm::ConsumesCollector& iCC);
  VertexTimeAlgorithmLegacy4D(edm::ParameterSet const& iConfig, edm::ConsumesCollector&& iCC) : VertexTimeAlgorithmLegacy4D(iConfig, iCC) {}
  ~VertexTimeAlgorithmLegacy4D() override = default;

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  void setEvent(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

  bool vertexTime(float& vtxTime, float& vtxTimeError, TransientVertex const& vtx) const override;
};

#endif
