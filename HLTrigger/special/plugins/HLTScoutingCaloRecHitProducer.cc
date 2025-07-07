#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingCaloRecHit.h"

class HLTScoutingCaloRecHitProducer : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingCaloRecHitProducer(const edm::ParameterSet&);
  ~HLTScoutingCaloRecHitProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& setup) const final;

  const edm::EDGetTokenT<reco::PFRecHitCollection> recoPFRecHitsToken_;
  const double minEnergy_;
  const int mantissaPrecision_;
};

HLTScoutingCaloRecHitProducer::HLTScoutingCaloRecHitProducer(const edm::ParameterSet& iConfig)
    : recoPFRecHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("pfRecHits"))),
      minEnergy_(iConfig.getParameter<double>("minEnergy")),
      mantissaPrecision_(iConfig.getParameter<int>("mantissaPrecision")) {
  produces<Run3ScoutingCaloRecHitCollection>();
}

void HLTScoutingCaloRecHitProducer::produce(edm::StreamID sid, edm::Event& iEvent, edm::EventSetup const& setup) const {
  auto const& recoPFRecHits = iEvent.get(recoPFRecHitsToken_);

  auto run3ScoutCaloRecHits = std::make_unique<Run3ScoutingCaloRecHitCollection>();
  run3ScoutCaloRecHits->reserve(recoPFRecHits.size());

  for (auto const& recoPFRecHit : recoPFRecHits) {
    if (recoPFRecHit.energy() < minEnergy_) {
      continue;
    }

    run3ScoutCaloRecHits->emplace_back(
        MiniFloatConverter::reduceMantissaToNbitsRounding(recoPFRecHit.energy(), mantissaPrecision_),
        recoPFRecHit.detId(),
        recoPFRecHit.flags());
  }

  iEvent.put(std::move(run3ScoutCaloRecHits));
}

void HLTScoutingCaloRecHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfRecHits", edm::InputTag("hltPFRecHits"));
  desc.add<double>("minEnergy", -1)->setComment("Minimum energy of the PFRecHit in GeV");
  desc.add<int>("mantissaPrecision", 10)->setComment("default of 10 corresponds to float16, change to 23 for float32");
  descriptions.add("hltScoutingCaloRecHitProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTScoutingCaloRecHitProducer);
