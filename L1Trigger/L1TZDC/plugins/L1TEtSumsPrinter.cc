#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1Trigger/interface/EtSum.h"

class L1TEtSumsPrinter : public edm::global::EDAnalyzer<> {
public:
  explicit L1TEtSumsPrinter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(edm::StreamID, edm::Event const&, edm::EventSetup const&) const override;

  edm::EDGetTokenT<l1t::EtSumBxCollection> srcToken_;
};

L1TEtSumsPrinter::L1TEtSumsPrinter(const edm::ParameterSet& iConfig)
    : srcToken_{consumes(iConfig.getParameter<edm::InputTag>("src"))} {}

void L1TEtSumsPrinter::analyze(edm::StreamID, edm::Event const& iEvent, edm::EventSetup const&) const {
  auto const& zdcEtSums = iEvent.get(srcToken_);
  auto const& moduleLabel = moduleDescription().moduleLabel();
  for (int ibx = zdcEtSums.getFirstBX(); ibx <= zdcEtSums.getLastBX(); ++ibx) {
    if (ibx != 0) {
      continue;
    }
    auto const size = zdcEtSums.size(ibx);
    for (uint idx = 0; idx < size; ++idx) {
      auto const& etSum = zdcEtSums.at(ibx, idx);

      auto const etSum_type = etSum.getType();

      if (etSum_type != l1t::EtSum::EtSumType::kMinBiasHFP0 and
          etSum_type != l1t::EtSum::EtSumType::kMinBiasHFM0 and
          etSum_type != l1t::EtSum::EtSumType::kMinBiasHFP1 and
          etSum_type != l1t::EtSum::EtSumType::kMinBiasHFM1) {
        continue;
      }

      edm::LogPrint("L1TEtSumsPrinter") << "[" << moduleLabel << "] EtSums[" << ibx << "][" << idx
                                        << "] (type, hwPt) = (" << etSum_type << ", " << etSum.hwPt() << ")";
    }
  }
}

void L1TEtSumsPrinter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("l1tEtSums"));
  descriptions.add("l1tEtSumsPrinter", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TEtSumsPrinter);
