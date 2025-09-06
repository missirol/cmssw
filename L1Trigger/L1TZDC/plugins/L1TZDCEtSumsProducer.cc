//
// L1TZDCEtSumsProducer
//  EDProducer to compute the ZDC EtSums from HCAL trigger primitives
//
// Original author: Chris McGinn
// Contact: christopher.mc.ginn@cern.ch or
//          cfmcginn on github for bugs/issues
//
#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

class L1TZDCEtSumsProducer : public edm::global::EDProducer<> {
public:
  explicit L1TZDCEtSumsProducer(edm::ParameterSet const& ps);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  edm::EDGetTokenT<HcalTrigPrimDigiCollection> const hcalTPDigisToken_;

  int const bxFirst_;
  int const bxLast_;

  static constexpr int kZDCAbsIEta = 42;
  static constexpr int kZDCEtSumsIPhi = 99;
  static constexpr int kZDCEtSumMaxValue = 1023;
};

L1TZDCEtSumsProducer::L1TZDCEtSumsProducer(edm::ParameterSet const& ps)
    : hcalTPDigisToken_{consumes(ps.getParameter<edm::InputTag>("hcalTPDigis"))},
      bxFirst_{ps.getParameter<int>("bxFirst")},
      bxLast_{ps.getParameter<int>("bxLast")} {
  produces<l1t::EtSumBxCollection>();
}

void L1TZDCEtSumsProducer::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& iSetup) const {
  auto outZDCEtSums = std::make_unique<l1t::EtSumBxCollection>(0, bxFirst_, bxLast_);

  auto const hcalTPs = iEvent.getHandle(hcalTPDigisToken_);

  if (not hcalTPs.isValid()) {
    edm::LogWarning("L1TZDCEtSumsProducer") << "Invalid handle to HcalTrigPrimDigiCollection collection"
                                            << ": returning empty l1t::EtSumBxCollection for ZDC EtSums !";
  } else if (bxFirst_ > bxLast_) {
    edm::LogWarning("L1TZDCEtSumsProducer")
        << "Invalid configuration parameters (bxFirst [" << bxFirst_ << "] > bxLast [" << bxLast_
        << "]): returning empty l1t::EtSumBxCollection for ZDC EtSums !";
  } else {
    unsigned int const nBXs = (bxLast_ - bxFirst_) + 1;

    std::vector<std::array<int, 2>> iEtSums{nBXs, {{0, 0}}};
    std::vector<std::array<bool, 2>> iEtSumsFillFlags{nBXs, {{false, false}}};

    for (auto const& hcalTp : *hcalTPs) {
      // iphi position 99 is used for the etSums
      auto const iphi = hcalTp.id().iphi();

      if (iphi == kZDCEtSumsIPhi) {
        continue;
      }

      // absIEta position 42 is used for the ZDC (-42 for ZDCM, +42 for ZDCP)
      auto const ieta = hcalTp.id().ieta();
      auto const absIEta = std::abs(ieta);

      if (absIEta != kZDCAbsIEta) {
        continue;
      }

      // index=0 (index=1) for ZDCM (ZDCP)
      auto const zdcIndex = (ieta < 0) ? 0 : 1;

      // get number of samples and number of presamples (nPresamples is bx=0)
      int const nSamples = hcalTp.size();
      int const nPresamples = hcalTp.presamples();

      for (auto iSample = 0; iSample < nSamples; ++iSample) {
        auto const ibx = iSample - nPresamples;
        if (ibx >= bxFirst_ and ibx <= bxLast_) {
          auto const& hcalTpSample = hcalTp.sample(iSample);
          auto const ietIn = hcalTpSample.raw() & kZDCEtSumMaxValue;
          auto const bxIndex = ibx - bxFirst_;
          iEtSums[bxIndex][zdcIndex] += ietIn;
          iEtSumsFillFlags[bxIndex][zdcIndex] = true;
        }
      }
    }

    for (unsigned int bxIndex = 0; bxIndex < iEtSums.size(); ++bxIndex) {
      int const bx = bxIndex + bxFirst_;
      for (unsigned int zdcIndex = 0; zdcIndex < 2; ++zdcIndex) {
        if (not iEtSumsFillFlags[bxIndex][zdcIndex]) {
          continue;
        }

        int const zdc_hwPt = std::min(iEtSums[bxIndex][zdcIndex], kZDCEtSumMaxValue);
        int const zdc_hwEta = (zdcIndex == 0) ? -1 : 1;

        auto const zdc_type = (zdcIndex == 0) ? l1t::EtSum::EtSumType::kZDCM : l1t::EtSum::EtSumType::kZDCP;

        l1t::EtSum zdc_etSum{};
        zdc_etSum.setHwPt(zdc_hwPt);
        zdc_etSum.setHwEta(zdc_hwEta);
        zdc_etSum.setHwPhi(0);
        zdc_etSum.setType(zdc_type);

        outZDCEtSums->push_back(bx, l1t::CaloTools::etSumP4Demux(zdc_etSum));
      }
    }
  }

  iEvent.put(std::move(outZDCEtSums));
}

void L1TZDCEtSumsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("hcalTPDigis", edm::InputTag("simHcalTriggerPrimitiveDigis"));
  desc.add<int>("bxFirst", -2);
  desc.add<int>("bxLast", 2);
  descriptions.add("l1tZDCEtSumsProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TZDCEtSumsProducer);
