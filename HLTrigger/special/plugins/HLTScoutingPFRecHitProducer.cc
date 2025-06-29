#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFRecHit.h"

class HLTScoutingPFRecHitProducer : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingPFRecHitProducer(const edm::ParameterSet&);
  ~HLTScoutingPFRecHitProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& setup) const final;

  const edm::EDGetTokenT<reco::PFRecHitCollection> recoPFRecHitsToken_;
  const double minEnergy_;
  const int mantissaPrecision_;
};

HLTScoutingPFRecHitProducer::HLTScoutingPFRecHitProducer(const edm::ParameterSet& iConfig)
    : recoPFRecHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
      minEnergy_(iConfig.getParameter<double>("minEnergy")),
      mantissaPrecision_(iConfig.getParameter<int>("mantissaPrecision")) {
  produces<Run3ScoutingPFRecHitCollection>();
}

void HLTScoutingPFRecHitProducer::produce(edm::StreamID sid, edm::Event& iEvent, edm::EventSetup const& setup) const {
  auto const& recoPFRecHits = iEvent.get(recoPFRecHitsToken_);

  auto run3ScoutPFRecHits = std::make_unique<Run3ScoutingPFRecHitCollection>();
  run3ScoutPFRecHits->reserve(recoPFRecHits.size());

  for (auto const& recoPFRecHit : recoPFRecHits) {
    if (recoPFRecHit.energy() < minEnergy_) {
      continue;
    }

    run3ScoutPFRecHits->emplace_back(
        recoPFRecHit.detId(),
        MiniFloatConverter::reduceMantissaToNbitsRounding(recoPFRecHit.energy(), mantissaPrecision_));
  }

  iEvent.put(std::move(run3ScoutPFRecHits));
}

void HLTScoutingPFRecHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltPFRecHits"));
  desc.add<double>("minEnergy", -1)->setComment("Minimum energy of the PFRecHit in GeV");
  desc.add<int>("mantissaPrecision", 10)->setComment("default of 10 corresponds to float16, change to 23 for float32");
  descriptions.add("hltScoutingPFRecHitProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTScoutingPFRecHitProducer);
