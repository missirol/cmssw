#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/JobConfigurationAlpakaRecord.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/alpaka/PFRecHitHBHEParamsAlpakaESData.h"

#include "TestAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TestPrintPFRecHitHBHEParamsOnDevice : public stream::EDProducer<> {
  public:
    TestPrintPFRecHitHBHEParamsOnDevice(edm::ParameterSet const& config) :
      esToken_{esConsumes(config.getParameter<edm::ESInputTag>("pfRecHitParams"))} {
      devicePutToken_ = produces("");
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      auto const& esData = iSetup.getData(esToken_);

      auto deviceProduct = std::make_unique<portabletest::TestDeviceCollection>(256, iEvent.queue());

      algo_.printPFRecHitHBHEParams(iEvent.queue(), esData);

      iEvent.put(devicePutToken_, std::move(deviceProduct));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::ESInputTag>("pfRecHitParams", edm::ESInputTag("pfRecHitHBHEParamsESProducer", ""));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::ESGetToken<PFRecHitHBHEParamsAlpakaESDataDevice, JobConfigurationAlpakaRecord> const esToken_;
    device::EDPutToken<portabletest::TestDeviceCollection> devicePutToken_;
    TestAlgo algo_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(TestPrintPFRecHitHBHEParamsOnDevice);
