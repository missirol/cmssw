#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "L1TriggerScouting/Utilities/interface/L1SCaloTowerRecoFixer.h"

class L1SCaloTowerReProcessor : public edm::global::EDProducer<> {
public:
  explicit L1SCaloTowerReProcessor(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static constexpr unsigned int kNBXPlus1 = 3565;

  edm::EDGetTokenT<l1ScoutingRun3::CaloTowerOrbitCollection> const src_;
  bool const fixBunchCrossing_;
};

L1SCaloTowerReProcessor::L1SCaloTowerReProcessor(const edm::ParameterSet& iPSet)
    : src_{consumes(iPSet.getParameter<edm::InputTag>("src"))},
      fixBunchCrossing_{iPSet.getParameter<bool>("fixBunchCrossing")} {
  produces<l1ScoutingRun3::CaloTowerOrbitCollection>();
}

void L1SCaloTowerReProcessor::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  auto const& inpCaloTowers = iEvent.get(src_);

  auto outCaloTowers = std::make_unique<l1ScoutingRun3::CaloTowerOrbitCollection>();
  std::vector<std::vector<l1ScoutingRun3::CaloTower>> outCaloTowersBuffer(kNBXPlus1);
  unsigned int nCaloTowersInOrbit = 0;

  for (auto const bx_old : inpCaloTowers.getFilledBxs()) {
    auto const& inpCaloTowersInThisBX = inpCaloTowers.bxIterator(bx_old);
    auto const nInpCaloTowersInThisBX = inpCaloTowersInThisBX.size();

    L1SCaloTowerRecoFixer fixer(iEvent.run(), iEvent.id().event(), bx_old);

    if (not(fixer.isFromMP70() or fixer.isFromMP71())) {
      throw cms::Exception("InvalidInput")
          << "Unexpected value for L1SCaloTowerRecoFixer::bx_correct() (its mod(9) must be equal to 7 or 8): "
          << fixer.bx_correct();
    }

    auto const bx_out = fixBunchCrossing_ ? fixer.bx_correct() : bx_old;

    auto& bufferThisBX = outCaloTowersBuffer[bx_out];
    bufferThisBX.reserve(nInpCaloTowersInThisBX);

    if (fixBunchCrossing_) {
      LogTrace("L1SCaloTowerReProcessor") << "[L1SCaloTowerReProcessor:" << moduleDescription().moduleLabel()
                                          << "] BX (old/input -> new/output) = " << bx_old << " -> " << bx_out;
    } else {
      LogTrace("L1SCaloTowerReProcessor")
          << "[L1SCaloTowerReProcessor:" << moduleDescription().moduleLabel() << "] BX = " << bx_old;
    }

    LogTrace("L1SCaloTowerReProcessor") << "[L1SCaloTowerReProcessor:" << moduleDescription().moduleLabel()
                                        << "]   L1-Scouting CaloTowers (old/input -> new/output)";

    for (auto idx = 0u; idx < nInpCaloTowersInThisBX; ++idx) {
      auto const& old_ct = inpCaloTowersInThisBX[idx];

      int16_t const new_hwEt = old_ct.hwEt();
      int16_t const new_erBits = old_ct.erBits();
      int16_t const new_miscBits = old_ct.miscBits();

      auto const correctCaloTowerHwEtaPhi = fixer.correctCaloTowerHwEtaAndHwPhi(old_ct.hwEta(), old_ct.hwPhi());
      int16_t const new_hwEta = correctCaloTowerHwEtaPhi.hwEta;
      int16_t const new_hwPhi = correctCaloTowerHwEtaPhi.hwPhi;

      LogTrace("L1SCaloTowerReProcessor")
          << "[L1SCaloTowerReProcessor:" << moduleDescription().moduleLabel() << "]     [" << idx
          << "] (hwEt=" << old_ct.hwEt() << ", hwEta=" << old_ct.hwEta() << ", hwPhi=" << old_ct.hwPhi()
          << ", erBits=" << old_ct.erBits() << ", miscBits=" << old_ct.miscBits() << ") -> (hwEt=" << new_hwEt
          << ", hwEta=" << new_hwEta << ", hwPhi=" << new_hwPhi << ", erBits=" << new_erBits
          << ", miscBits=" << new_miscBits << ")";

      bufferThisBX.emplace_back(new_hwEt, new_erBits, new_miscBits, new_hwEta, new_hwPhi);
      ++nCaloTowersInOrbit;
    }
  }

  outCaloTowers->fillAndClear(outCaloTowersBuffer, nCaloTowersInOrbit);

  iEvent.put(std::move(outCaloTowers));
}

void L1SCaloTowerReProcessor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment(
      "Input collection of CaloTowers (type: l1ScoutingRun3::CaloTowerOrbitCollection)");
  desc.add<bool>("fixBunchCrossing", false)
      ->setComment("If the bunch crossing value is incorrect, use the correct one");
  descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1SCaloTowerReProcessor);
