#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalObjectMapRecord.h"

#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

class L1TAcceptAnalyzer : public edm::global::EDAnalyzer<> {
public:
  explicit L1TAcceptAnalyzer(const edm::ParameterSet&);
  ~L1TAcceptAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void analyze(edm::StreamID, edm::Event const&, const edm::EventSetup&) const override;

private:
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> const globalAlgBlkBx_token_;
  edm::EDGetTokenT<GlobalObjectMapRecord> const globalObjMapRecord_token_;
};

L1TAcceptAnalyzer::L1TAcceptAnalyzer(const edm::ParameterSet& iPSet)
  : globalAlgBlkBx_token_(consumes<GlobalAlgBlkBxCollection>(iPSet.getParameter<edm::InputTag>("globalAlgoBlocks"))),
    globalObjMapRecord_token_(consumes<GlobalObjectMapRecord>(iPSet.getParameter<edm::InputTag>("globalObjectMapRecord"))) {}

void L1TAcceptAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("globalAlgoBlocks", edm::InputTag("hltGtStage2Digis"));
  desc.add<edm::InputTag>("globalObjectMapRecord", edm::InputTag("hltGtStage2ObjectMap"));
  descriptions.add("l1tAcceptAnalyzer", desc);
}

void L1TAcceptAnalyzer::analyze(edm::StreamID, edm::Event const& iEvent, edm::EventSetup const& iSetup) const {

  // get handle to unpacked GT
  auto const& uGtAlgoBlocks = iEvent.get(globalAlgBlkBx_token_);

  // check size (BX 0)
  if (uGtAlgoBlocks.isEmpty(0)) {
    edm::LogWarning("L1TAcceptAnalyzer") << " Warning: GlobalAlgBlkBxCollection is empty for BX=0.";
    return;
  }

  // get handle to object maps from emulator (one object map per algorithm)
  auto const& gtObjectMapRecord = iEvent.get(globalObjMapRecord_token_);
  auto const& objMaps = gtObjectMapRecord.gtObjectMap();

  for (size_t imap = 0; imap < objMaps.size(); imap++) {
    auto const bit = objMaps[imap].algoBitNumber();
    auto const initDecision = uGtAlgoBlocks.at(0, 0).getAlgoDecisionInitial(bit);
    auto const emulDecision = objMaps[imap].algoGtlResult();

    if (emulDecision != initDecision) {
      edm::LogPrint("L1TAcceptAnalyzer") << iEvent.id() << " -- L1T decision (emulated vs. unpacked initial) is not the same: "
                                         << objMaps[imap].algoName() << " (bit = " << bit << "), emulated = "
                                         << emulDecision << ", unpacked (initial) = " << initDecision;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TAcceptAnalyzer);
