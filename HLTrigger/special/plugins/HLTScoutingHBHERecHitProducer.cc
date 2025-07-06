#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingHBHERecHit.h"

class HLTScoutingHBHERecHitProducer : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingHBHERecHitProducer(const edm::ParameterSet&);
  ~HLTScoutingHBHERecHitProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const final;

  const edm::EDGetTokenT<reco::PFRecHitCollection> recoPFRecHitsToken_;
  const double minEnergy_;
  const int mantissaPrecision_;
};

HLTScoutingHBHERecHitProducer::HLTScoutingHBHERecHitProducer(const edm::ParameterSet& iConfig)
    : recoPFRecHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("pfRecHits"))),
      minEnergy_(iConfig.getParameter<double>("minEnergy")),
      mantissaPrecision_(iConfig.getParameter<int>("mantissaPrecision")) {
  produces<Run3ScoutingHBHERecHitCollection>();
}

void HLTScoutingHBHERecHitProducer::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  auto const& recoPFRecHits = iEvent.get(recoPFRecHitsToken_);

  auto run3ScoutHBHERecHits = std::make_unique<Run3ScoutingHBHERecHitCollection>();
  run3ScoutHBHERecHits->reserve(recoPFRecHits.size());

  for (auto const& recoPFRecHit : recoPFRecHits) {
    if (recoPFRecHit.energy() < minEnergy_) {
      continue;
    }

    HcalDetId const hcalDetId{recoPFRecHit.detId()};
    run3ScoutHBHERecHits->emplace_back(
        MiniFloatConverter::reduceMantissaToNbitsRounding(recoPFRecHit.energy(), mantissaPrecision_),
        hcalDetId.ieta(),
        hcalDetId.iphi(),
        hcalDetId.depth());
  }

  iEvent.put(std::move(run3ScoutHBHERecHits));
}

void HLTScoutingHBHERecHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfRecHits", edm::InputTag("hltPFRecHits"));
  desc.add<double>("minEnergy", -1)->setComment("Minimum energy of the PFRecHit in GeV");
  desc.add<int>("mantissaPrecision", 10)->setComment("default of 10 corresponds to float16, change to 23 for float32");
  descriptions.add("hltScoutingHBHERecHitProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTScoutingHBHERecHitProducer);
