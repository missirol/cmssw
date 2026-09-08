#include <memory>
#include <utility>
#include <vector>

#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloJet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "L1TriggerScouting/OnlineProcessing/interface/L1ScoutingCaloJetClusterizer.h"

class L1ScoutingCaloJetOrbitProducer : public edm::global::EDProducer<> {
public:
  explicit L1ScoutingCaloJetOrbitProducer(const edm::ParameterSet&);
  ~L1ScoutingCaloJetOrbitProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  // Number of BXs per orbit plus one
  static constexpr unsigned int kNBXPlus1 = 3565;

  edm::EDGetTokenT<l1ScoutingRun3::CaloTowerOrbitCollection> const src_;
  bool const produceSortedCaloTowers_;
  L1ScoutingCaloJetClusterizer const l1sCaloJetClusterizer_;
};

L1ScoutingCaloJetOrbitProducer::L1ScoutingCaloJetOrbitProducer(const edm::ParameterSet& iPSet)
    : src_(consumes(iPSet.getParameter<edm::InputTag>("src"))),
      produceSortedCaloTowers_{iPSet.getParameter<bool>("produceSortedCaloTowers")},
      l1sCaloJetClusterizer_{iPSet} {
  produces<l1ScoutingRun3::CaloJetOrbitCollection>("CaloJet").setBranchAlias("CaloJetOrbitCollection");
  if (produceSortedCaloTowers_) {
    produces<l1ScoutingRun3::CaloTowerOrbitCollection>("SortedCaloTowers");
  }
}

void L1ScoutingCaloJetOrbitProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  // Input collection of CaloTowers
  auto const& caloTowerCollection{iEvent.get(src_)};

  // Output containers for CaloJets
  auto caloJetCollection{std::make_unique<l1ScoutingRun3::CaloJetOrbitCollection>()};
  std::vector<std::vector<l1ScoutingRun3::CaloJet>> caloJetBuffer(kNBXPlus1);
  auto nCaloJet{0u};

  // Output containers for sorted CaloTowers (used only if "produceSortedCaloTowers == True")
  auto sortedCaloTowerCollection{std::make_unique<l1ScoutingRun3::CaloTowerOrbitCollection>()};
  std::vector<std::vector<l1ScoutingRun3::CaloTower>> sortedCaloTowerBuffer(kNBXPlus1);
  auto nSortedCaloTower{0u};

  auto const& moduleName{moduleDescription().moduleName()};
  auto const& moduleLabel{moduleDescription().moduleLabel()};

  for (auto const bx : caloTowerCollection.getFilledBxs()) {
    LogTrace(moduleName) << "[" << moduleName << ":" << moduleLabel << "] BX = " << bx;
    LogTrace(moduleName) << "[" << moduleName << ":" << moduleLabel
                         << "]   Inputs (l1ScoutingRun3::CaloTower and fastjet::PseudoJet)";

    l1sCaloJetClusterizer_.run(
        caloTowerCollection.bxIterator(bx), caloJetBuffer[bx], sortedCaloTowerBuffer[bx], produceSortedCaloTowers_);

    nCaloJet += caloJetBuffer[bx].size();
    nSortedCaloTower += sortedCaloTowerBuffer[bx].size();
  }

  // Fill orbit collection(s) with output product(s)
  caloJetCollection->fillAndClear(caloJetBuffer, nCaloJet);
  iEvent.put(std::move(caloJetCollection), "CaloJet");

  if (produceSortedCaloTowers_) {
    sortedCaloTowerCollection->fillAndClear(sortedCaloTowerBuffer, nSortedCaloTower);
    iEvent.put(std::move(sortedCaloTowerCollection), "SortedCaloTowers");
  }

  LogTrace(moduleName) << "[" << moduleName << ":" << moduleLabel << "] === End of produce() method ==";
}

void L1ScoutingCaloJetOrbitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src")->setComment(
      "Input collection of CaloTowers (type: l1ScoutingRun3::CaloTowerOrbitCollection)");

  desc.add<bool>("produceSortedCaloTowers", false)
      ->setComment("Output a copy of the l1ScoutingRun3::CaloTowerOrbitCollection in \"src\" with a custom sorting");

  L1ScoutingCaloJetClusterizer::fillDescription(desc);

  descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ScoutingCaloJetOrbitProducer);
