#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1Trigger/interface/EtSum.h"

class L1TZDCEtSumsPrinter : public edm::global::EDAnalyzer<> {
public:
  explicit L1TZDCEtSumsPrinter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(edm::StreamID, edm::Event const&, edm::EventSetup const&) const override;

  edm::EDGetTokenT<l1t::EtSumBxCollection> srcToken_;
};

L1TZDCEtSumsPrinter::L1TZDCEtSumsPrinter(const edm::ParameterSet& iConfig)
    : srcToken_{consumes(iConfig.getParameter<edm::InputTag>("src"))} {}

void L1TZDCEtSumsPrinter::analyze(edm::StreamID, edm::Event const& iEvent, edm::EventSetup const&) const {
  auto const& zdcEtSums = iEvent.get(srcToken_);
  auto const& moduleLabel = moduleDescription().moduleLabel();
  for (int ibx = zdcEtSums.getFirstBX(); ibx <= zdcEtSums.getLastBX(); ++ibx) {
    auto const size = zdcEtSums.size(ibx);
    for (uint idx = 0; idx < size; ++idx) {
      auto const& etSum = zdcEtSums.at(ibx, idx);
      edm::LogPrint("L1TZDCEtSumsPrinter") << "[" << moduleLabel << "] zdcEtSums[" << ibx << "][" << idx
                                           << "] (type, hwPt) = (" << etSum.getType() << ", " << etSum.hwPt() << ")";
    }
  }
}

void L1TZDCEtSumsPrinter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("l1tZDCEtSums"));
  descriptions.add("l1tZDCEtSumsPrinter", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TZDCEtSumsPrinter);
