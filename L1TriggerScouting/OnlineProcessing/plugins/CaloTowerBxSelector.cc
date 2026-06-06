#include <memory>
#include <utility>
#include <vector>

#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

class CaloTowerBxSelector : public edm::global::EDProducer<> {
public:
  explicit CaloTowerBxSelector(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  edm::EDGetTokenT<l1ScoutingRun3::CaloTowerOrbitCollection> const towersToken_;
  unsigned int const minNTower_;
};

CaloTowerBxSelector::CaloTowerBxSelector(const edm::ParameterSet& iPSet)
    : towersToken_(consumes(iPSet.getParameter<edm::InputTag>("towersTag"))),
      minNTower_(iPSet.getParameter<unsigned int>("minNTower")) {
  produces<std::vector<unsigned int>>("SelBx").setBranchAlias("CaloTowerSelectedBx");
}

void CaloTowerBxSelector::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  auto const& caloTowersInOrbit = iEvent.get(towersToken_);

  auto towerBx = std::make_unique<std::vector<unsigned int>>();

  // loop over valid BXs with CaloTowers
  for (auto const bx : caloTowersInOrbit.getFilledBxs()) {
    if (minNTower_ > 0) {
      auto const& towers = caloTowersInOrbit.bxIterator(bx);

      // skip BX if the number of towers is below the minimum
      if (towers.size() < minNTower_) {
        continue;
      }
    }

    towerBx->emplace_back(bx);
  }  // end orbit loop

  iEvent.put(std::move(towerBx), "SelBx");
}

void CaloTowerBxSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("towersTag")
      ->setComment("Input collection of CaloTowers (type: l1ScoutingRun3::CaloTowerOrbitCollection)");
  desc.add<unsigned int>("minNTower", 0)->setComment("Min number of towers");
  descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CaloTowerBxSelector);
