#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmBase.h"

void VertexTimeAlgorithmBase::fill_vertex_times(std::vector<TransientVertex>& pvs) const {
  for (unsigned int idx = 0; idx < pvs.size(); ++idx) {
    auto const& vtx = pvs[idx];

    if (not vtx.isValid()) {
      continue;
    }

    auto vtxTime(0.f), vtxTimeError(-1.f);
    if (not vertexTime(vtxTime, vtxTimeError, vtx)) {
      continue;
    }

    auto err = vtx.positionError().matrix4D();
    err(3, 3) = vtxTimeError * vtxTimeError;
    auto vtx_with_time = TransientVertex(vtx.position(), vtxTime, err, vtx.originalTracks(), vtx.totalChiSquared(), vtx.degreesOfFreedom());
    vtx_with_time.weightMap(vtx.weightMap());
    pvs[idx] = vtx_with_time;
  }
}
