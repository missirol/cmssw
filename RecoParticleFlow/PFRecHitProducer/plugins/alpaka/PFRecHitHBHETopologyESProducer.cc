#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/JobConfigurationAlpakaRecord.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/alpaka/PFRecHitHBHETopologyAlpakaESData.h"

//#include "RecoParticleFlow/PFClusterProducer/interface/PFHCALDenseIdNavigatorCore.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PFRecHitHBHETopologyESProducer : public ESProducer {
  public:
    PFRecHitHBHETopologyESProducer(edm::ParameterSet const& iConfig) {
      setWhatProduced(this);
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<std::string>("appendToDataLabel", "");
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<PFRecHitHBHETopologyAlpakaESDataHost> produce(JobConfigurationAlpakaRecord const& iRecord) {
      uint32_t const productSize = 1; //!!

      auto product = std::make_unique<PFRecHitHBHETopologyAlpakaESDataHost>(productSize, cms::alpakatools::host());

      fillPFRecHitHBHETopologyAlpakaESDataHost(product->view());

      return product;
    }

  private:
    void fillPFRecHitHBHETopologyAlpakaESDataHost(PFRecHitHBHETopologyAlpakaESDataHost::View view) const;
  };

  void PFRecHitHBHETopologyESProducer::fillPFRecHitHBHETopologyAlpakaESDataHost(PFRecHitHBHETopologyAlpakaESDataHost::View view) const {
    // FILL
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(PFRecHitHBHETopologyESProducer);
