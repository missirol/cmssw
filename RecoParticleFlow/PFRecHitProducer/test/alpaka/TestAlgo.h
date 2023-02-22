#ifndef RecoParticleFlow_PFRecHitProducer_plugins_alpaka_TestAlgo_h
#define RecoParticleFlow_PFRecHitProducer_plugins_alpaka_TestAlgo_h

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/alpaka/PFRecHitHBHEParamsAlpakaESData.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TestAlgo {
  public:
    void printPFRecHitHBHEParams(Queue& queue, PFRecHitHBHEParamsAlpakaESDataDevice const&) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoParticleFlow_PFRecHitProducer_plugins_alpaka_TestAlgo_h
