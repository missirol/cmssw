#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

class L1ScoutingCaloTowerPhysicalValueMapProducer : public edm::global::EDProducer<> {
public:
  L1ScoutingCaloTowerPhysicalValueMapProducer(edm::ParameterSet const&);
  ~L1ScoutingCaloTowerPhysicalValueMapProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  void putValueMap(edm::Event&,
                   edm::Handle<l1ScoutingRun3::CaloTowerOrbitCollection> const&,
                   std::vector<float> const&,
                   std::string const&) const;

  edm::EDGetTokenT<l1ScoutingRun3::CaloTowerOrbitCollection> const src_;
};

L1ScoutingCaloTowerPhysicalValueMapProducer::L1ScoutingCaloTowerPhysicalValueMapProducer(edm::ParameterSet const& params)
    : src_(consumes(params.getParameter<edm::InputTag>("src"))) {
  produces<edm::ValueMap<float>>("fEt");
  produces<edm::ValueMap<float>>("fEta");
  produces<edm::ValueMap<float>>("fPhi");
}

void L1ScoutingCaloTowerPhysicalValueMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  descriptions.addWithDefaultLabel(desc);
}

void L1ScoutingCaloTowerPhysicalValueMapProducer::produce(edm::StreamID,
                                                          edm::Event& iEvent,
                                                          edm::EventSetup const&) const {
  auto const src_h = iEvent.getHandle(src_);
  auto const& src = *src_h;
  auto const nobjs = src.size();

  std::vector<float> outv_fEt{};
  outv_fEt.reserve(nobjs);

  std::vector<float> outv_fEta{};
  outv_fEta.reserve(nobjs);

  std::vector<float> outv_fPhi{};
  outv_fPhi.reserve(nobjs);

  for (auto iobj = 0; iobj < nobjs; ++iobj) {
    auto const& obj = src[iobj];
    // conversion from (hwEt, hwEta, hwPhi) to (fEt, fEta, fPhi) taken from l1t::CaloTools::p4MP
    outv_fEt.emplace_back(obj.hwEt() * 0.5 + 1E-6);
    outv_fEta.emplace_back(l1t::CaloTools::towerEta(obj.hwEta()));
    outv_fPhi.emplace_back(l1t::CaloTools::towerPhi(obj.hwEta(), obj.hwPhi()));
  }

  putValueMap(iEvent, src_h, outv_fEt, "fEt");
  putValueMap(iEvent, src_h, outv_fEta, "fEta");
  putValueMap(iEvent, src_h, outv_fPhi, "fPhi");
}

void L1ScoutingCaloTowerPhysicalValueMapProducer::putValueMap(
    edm::Event& iEvent,
    edm::Handle<l1ScoutingRun3::CaloTowerOrbitCollection> const& handle,
    std::vector<float> const& values,
    std::string const& label) const {
  auto valuemap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler filler(*valuemap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(valuemap), label);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ScoutingCaloTowerPhysicalValueMapProducer);
