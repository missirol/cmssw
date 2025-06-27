#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/Scouting/interface/Run3ScoutingCaloTower2.h"

class HLTScoutingCaloTower2Producer : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingCaloTower2Producer(const edm::ParameterSet&);
  ~HLTScoutingCaloTower2Producer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& setup) const final;

  const edm::EDGetTokenT<CaloTowerCollection> recoCaloTowersToken_;
  const double minEnergy_;
  const int mantissaPrecision_;
};

HLTScoutingCaloTower2Producer::HLTScoutingCaloTower2Producer(const edm::ParameterSet& iConfig)
    : recoCaloTowersToken_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
      minEnergy_(iConfig.getParameter<double>("minEnergy")),
      mantissaPrecision_(iConfig.getParameter<int>("mantissaPrecision")) {
  produces<Run3ScoutingCaloTower2Collection>();
}

void HLTScoutingCaloTower2Producer::produce(edm::StreamID sid, edm::Event& iEvent, edm::EventSetup const& setup) const {

  auto const& recoCaloTowers = iEvent.get(recoCaloTowersToken_);

  auto run3ScoutCaloTowers = std::make_unique<Run3ScoutingCaloTower2Collection>();
  run3ScoutCaloTowers->reserve(recoCaloTowers.size());

  for (auto const& recoCaloTower : recoCaloTowers) {

    if (recoCaloTower.energy() < minEnergy_) {
      continue;
    }

    run3ScoutCaloTowers->emplace_back(
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.emEnergy(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.hadEnergy(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.outerEnergy(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.p4().eta(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.p4().phi(), mantissaPrecision_)
    );
  }

  iEvent.put(std::move(run3ScoutCaloTowers));
}

void HLTScoutingCaloTower2Producer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltTowerMakerForAll"));
  desc.add<double>("minEnergy", -1)->setComment("Minimum energy of the CaloTower in GeV");
  desc.add<int>("mantissaPrecision", 10)->setComment("default of 10 corresponds to float16, change to 23 for float32");
  descriptions.add("hltScoutingCaloTower2Producer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTScoutingCaloTower2Producer);
