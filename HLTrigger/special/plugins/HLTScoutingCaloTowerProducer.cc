#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/Scouting/interface/Run3ScoutingCaloTower.h"

class HLTScoutingCaloTowerProducer : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingCaloTowerProducer(const edm::ParameterSet&);
  ~HLTScoutingCaloTowerProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& setup) const final;

  const edm::EDGetTokenT<CaloTowerCollection> recoCaloTowersToken_;
  const int mantissaPrecision_;
};

HLTScoutingCaloTowerProducer::HLTScoutingCaloTowerProducer(const edm::ParameterSet& iConfig)
    : recoCaloTowersToken_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
      mantissaPrecision_(iConfig.getParameter<int>("mantissaPrecision")) {
  produces<Run3ScoutingCaloTowerCollection>();
}

void HLTScoutingCaloTowerProducer::produce(edm::StreamID sid, edm::Event& iEvent, edm::EventSetup const& setup) const {

  auto const& recoCaloTowers = iEvent.get(recoCaloTowersToken_);

  auto run3ScoutCaloTowers = std::make_unique<Run3ScoutingCaloTowerCollection>();
  run3ScoutCaloTowers->reserve(recoCaloTowers.size());

  for (auto const& recoCaloTower : recoCaloTowers) {

    run3ScoutCaloTowers->emplace_back(
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.p4().pt(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.p4().eta(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.p4().phi(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.p4().mass(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.emEnergy(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.hadEnergy(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.outerEnergy(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.ecalTime(), mantissaPrecision_),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloTower.hcalTime(), mantissaPrecision_),
      recoCaloTower.ieta(),
      recoCaloTower.iphi(),
      recoCaloTower.numCrystals(),
      recoCaloTower.constituents().size(),
      recoCaloTower.towerStatusWord()
    );
  }

  iEvent.put(std::move(run3ScoutCaloTowers));
}

void HLTScoutingCaloTowerProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltTowerMakerForAll"));
  desc.add<int>("mantissaPrecision", 10)->setComment("default float16, change to 23 for float32");
  descriptions.add("hltScoutingCaloTowerProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTScoutingCaloTowerProducer);
