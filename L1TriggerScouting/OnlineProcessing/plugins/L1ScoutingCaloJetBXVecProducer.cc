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
#include "L1TriggerScouting/Utilities/interface/L1ScoutingBXVectors.h"

class L1ScoutingCaloJetBXVecProducer : public edm::global::EDProducer<> {
public:
  explicit L1ScoutingCaloJetBXVecProducer(const edm::ParameterSet&);
  ~L1ScoutingCaloJetBXVecProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  edm::EDGetTokenT<l1ScoutingRun3::CaloTowerBxCollection> const src_;
  bool const produceSortedCaloTowers_;
  L1ScoutingCaloJetClusterizer const l1sCaloJetClusterizer_;
};

L1ScoutingCaloJetBXVecProducer::L1ScoutingCaloJetBXVecProducer(const edm::ParameterSet& iPSet)
    : src_(consumes(iPSet.getParameter<edm::InputTag>("src"))),
      produceSortedCaloTowers_{iPSet.getParameter<bool>("produceSortedCaloTowers")},
      l1sCaloJetClusterizer_{iPSet} {
  produces<l1ScoutingRun3::CaloJetBxCollection>("CaloJets");
  if (produceSortedCaloTowers_) {
    produces<l1ScoutingRun3::CaloTowerBxCollection>("SortedCaloTowers");
  }
}

void L1ScoutingCaloJetBXVecProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  auto const& input{iEvent.get(src_)};

  auto const bxMin{input.getFirstBX()};
  auto const bxMax{input.getLastBX()};

  auto outputCaloJets{std::make_unique<l1ScoutingRun3::CaloJetBxCollection>(0, bxMin, bxMax)};

  // Output containers for sorted CaloTowers (used only if "produceSortedCaloTowers == True")
  auto outputCaloTowers{std::make_unique<l1ScoutingRun3::CaloTowerBxCollection>(0, bxMin, bxMax)};

  auto const& moduleName{moduleDescription().moduleName()};
  auto const& moduleLabel{moduleDescription().moduleLabel()};

  for (auto bx = bxMin; bx <= bxMax; ++bx) {
    LogTrace(moduleName) << "[" << moduleName << ":" << moduleLabel << "] BX = " << bx;
    LogTrace(moduleName) << "[" << moduleName << ":" << moduleLabel
                         << "]   Inputs (l1ScoutingRun3::CaloTower and fastjet::PseudoJet)";

    std::vector<l1ScoutingRun3::CaloJet> bxCaloJetBuffer{};
    std::vector<l1ScoutingRun3::CaloTower> bxCaloTowerBuffer{};

    l1sCaloJetClusterizer_.run(input.bxIterator(bx), bxCaloJetBuffer, bxCaloTowerBuffer, produceSortedCaloTowers_);

    for (auto const& caloJet : bxCaloJetBuffer) {
      outputCaloJets->push_back(bx, caloJet);
    }

    if (produceSortedCaloTowers_) {
      for (auto const& caloTower : bxCaloTowerBuffer) {
        outputCaloTowers->push_back(bx, caloTower);
      }
    }
  }

  iEvent.put(std::move(outputCaloJets), "CaloJets");

  if (produceSortedCaloTowers_) {
    iEvent.put(std::move(outputCaloTowers), "SortedCaloTowers");
  }

  LogTrace(moduleName) << "[" << moduleName << ":" << moduleLabel << "] === End of produce() method ==";
}

void L1ScoutingCaloJetBXVecProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src")->setComment(
      "Input collection of CaloTowers (type: l1ScoutingRun3::CaloTowerBxCollection)");

  desc.add<bool>("produceSortedCaloTowers", false)
      ->setComment("Output a copy of the l1ScoutingRun3::CaloTowerBxCollection in \"src\" with a custom sorting");

  L1ScoutingCaloJetClusterizer::fillDescription(desc);

  descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ScoutingCaloJetBXVecProducer);
