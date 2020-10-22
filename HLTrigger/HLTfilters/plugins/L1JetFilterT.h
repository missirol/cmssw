#ifndef HLTrigger_HLTfilters_L1JetFilterT_h
#define HLTrigger_HLTfilters_L1JetFilterT_h

#include <vector>
#include <cmath>

#include "CondFormats/HLTObjects/interface/L1TObjScalingConstants.h"
#include "CondFormats/DataRecord/interface/L1TObjScalingRcd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"

template<class T>
class L1JetFilterT : public HLTFilter {
public:
  explicit L1JetFilterT(const edm::ParameterSet&);
  ~L1JetFilterT() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs& filterproduct) const override;

private:
  edm::InputTag l1tJetTag_;
  edm::EDGetTokenT<std::vector<T>> l1tJetToken_;
  edm::ESInputTag l1tJetScalingTag_;
  edm::ESGetToken<L1TObjScalingConstants, L1TObjScalingRcd> l1tJetScalingToken_;

  double const minPt_;
  double const minEta_;
  double const maxEta_;
  int const minN_;

  double offlineJetPt(double const pt, double const eta, L1TObjScalingConstants const& scalingConstants) const;
};

template<class T>
L1JetFilterT<T>::L1JetFilterT(const edm::ParameterSet& iConfig)
    : HLTFilter(iConfig),
      l1tJetTag_(iConfig.getParameter<edm::InputTag>("inputTag")),
      l1tJetToken_(consumes<std::vector<T>>(l1tJetTag_)),
      l1tJetScalingTag_(iConfig.getParameter<edm::ESInputTag>("esScalingTag")),
      l1tJetScalingToken_(esConsumes<L1TObjScalingConstants, L1TObjScalingRcd>(l1tJetScalingTag_)),
      minPt_(iConfig.getParameter<double>("MinPt")),
      minEta_(iConfig.getParameter<double>("MinEta")),
      maxEta_(iConfig.getParameter<double>("MaxEta")),
      minN_(iConfig.getParameter<int>("MinN")) {}

template<class T>
L1JetFilterT<T>::~L1JetFilterT() = default;

template<class T>
void L1JetFilterT<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<edm::InputTag>("inputTag", edm::InputTag("ak4PFL1PuppiCorrected"));
  desc.add<edm::ESInputTag>("esScalingTag", edm::ESInputTag("L1TScalingESSource", "L1PFJetScaling"));
  desc.add<double>("MinPt", -1.0);
  desc.add<double>("MinEta", -5.0);
  desc.add<double>("MaxEta", 5.0);
  desc.add<int>("MinN", 1);
  descriptions.add(defaultModuleLabel<L1JetFilterT<T>>(), desc);
}

template<class T>
bool L1JetFilterT<T>::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs& filterproduct) const {
  // All HLT filters must create and fill an HLT filter object,
  // recording any reconstructed physics objects satisfying (or not)
  // this HLT filter, and place it in the Event.

  // The filter object
  if (saveTags()) {
    filterproduct.addCollectionTag(l1tJetTag_);
  }

  auto const& l1tJets = iEvent.getHandle(l1tJetToken_);
  auto const& l1tScalingConstants = iSetup.getData(l1tJetScalingToken_);

  int nJet(0);
  for (auto iJet = l1tJets->begin(); iJet != l1tJets->end(); ++iJet) {
    if (offlineJetPt(iJet->pt(), iJet->eta(), l1tScalingConstants) >= minPt_ && iJet->eta() <= maxEta_ && iJet->eta() >= minEta_) {
      ++nJet;
      edm::Ref<std::vector<T>> ref(l1tJets, std::distance(l1tJets->begin(), iJet));
      filterproduct.addObject(trigger::TriggerObjectType::TriggerL1PFJet, ref);
    }
  }

  // return with final filter decision
  return nJet >= minN_;
}

template<class T>
double L1JetFilterT<T>::offlineJetPt(double const pt, double const eta, const L1TObjScalingConstants& scalingConstants) const {
  uint const scIdx = std::abs(eta) < 1.5 ? 0 : (std::abs(eta) < 2.4 ? 1 : 2);

  if (scIdx >= scalingConstants.m_constants.size()) {
    throw cms::Exception("Input") << "out-of-range index for L1TObjScalingConstants vector (size=" << scalingConstants.m_constants.size() << ") [jet: pt=" << pt << ", eta="<< eta << "]";
  }

  return scalingConstants.m_constants.at(scIdx).m_constant + pt * scalingConstants.m_constants.at(scIdx).m_linear + pt * pt * scalingConstants.m_constants.at(scIdx).m_quadratic;
}

#endif // HLTrigger_HLTfilters_L1JetFilterT_h
