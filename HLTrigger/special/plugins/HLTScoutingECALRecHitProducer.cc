#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingEBRecHit.h"
#include "DataFormats/Scouting/interface/Run3ScoutingEERecHit.h"

class HLTScoutingECALRecHitProducer : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingECALRecHitProducer(const edm::ParameterSet&);
  ~HLTScoutingECALRecHitProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const final;

  const edm::EDGetTokenT<reco::PFRecHitCollection> recoPFRecHitsToken_;
  const double minEnergyEB_;
  const double minEnergyEE_;
  const int mantissaPrecision_;
};

HLTScoutingECALRecHitProducer::HLTScoutingECALRecHitProducer(const edm::ParameterSet& iConfig)
    : recoPFRecHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("pfRecHits"))),
      minEnergyEB_(iConfig.getParameter<double>("minEnergyEB")),
      minEnergyEE_(iConfig.getParameter<double>("minEnergyEE")),
      mantissaPrecision_(iConfig.getParameter<int>("mantissaPrecision")) {
  produces<Run3ScoutingEBRecHitCollection>("EB");
  produces<Run3ScoutingEERecHitCollection>("EE");
}

void HLTScoutingECALRecHitProducer::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  auto const& recoPFRecHits = iEvent.get(recoPFRecHitsToken_);

  auto run3ScoutEBRecHits = std::make_unique<Run3ScoutingEBRecHitCollection>();
  run3ScoutEBRecHits->reserve(recoPFRecHits.size());

  auto run3ScoutEERecHits = std::make_unique<Run3ScoutingEERecHitCollection>();
  run3ScoutEERecHits->reserve(recoPFRecHits.size());

  for (auto const& recoPFRecHit : recoPFRecHits) {
    if (recoPFRecHit.layer() == PFLayer::ECAL_BARREL) {
      if (recoPFRecHit.energy() < minEnergyEB_) {
        continue;
      }

      EBDetId const ebDetId{recoPFRecHit.detId()};
      run3ScoutEBRecHits->emplace_back(
          MiniFloatConverter::reduceMantissaToNbitsRounding(recoPFRecHit.energy(), mantissaPrecision_),
          ebDetId.ieta(),
          ebDetId.iphi());
    } else if (recoPFRecHit.layer() == PFLayer::ECAL_ENDCAP) {
      if (recoPFRecHit.energy() < minEnergyEE_) {
        continue;
      }

      EEDetId const eeDetId{recoPFRecHit.detId()};
      run3ScoutEERecHits->emplace_back(
          MiniFloatConverter::reduceMantissaToNbitsRounding(recoPFRecHit.energy(), mantissaPrecision_),
          eeDetId.ix(),
          eeDetId.iy(),
          eeDetId.positiveZ());
    } else {
      edm::LogWarning("HLTScoutingECALRecHitProducer")
          << "Skipping PFRecHit because of unexpected PFLayer value (" << recoPFRecHit.layer() << ").";
    }
  }

  iEvent.put(std::move(run3ScoutEBRecHits), "EB");
  iEvent.put(std::move(run3ScoutEERecHits), "EE");
}

void HLTScoutingECALRecHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfRecHits", edm::InputTag("hltPFRecHits"));
  desc.add<double>("minEnergyEB", -1)->setComment("Minimum energy of the EcalBarrel PFRecHit in GeV");
  desc.add<double>("minEnergyEE", -1)->setComment("Minimum energy of the EcalEndcap PFRecHit in GeV");
  desc.add<int>("mantissaPrecision", 10)->setComment("default of 10 corresponds to float16, change to 23 for float32");
  descriptions.add("hltScoutingECALRecHitProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTScoutingECALRecHitProducer);
