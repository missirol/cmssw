#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/Scouting/interface/Run3ScoutingCaloRecHit.h"

template<typename T>
class HLTScoutingCaloRecHitProducerT : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingCaloRecHitProducerT(const edm::ParameterSet&);
  ~HLTScoutingCaloRecHitProducerT() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& setup) const final;

  const edm::EDGetTokenT<T> recoCaloRecHitsToken_;
  const double minEnergy_;
  const int mantissaPrecision_;
};

template<typename T>
HLTScoutingCaloRecHitProducerT<T>::HLTScoutingCaloRecHitProducerT(const edm::ParameterSet& iConfig)
    : recoCaloRecHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
      minEnergy_(iConfig.getParameter<double>("minEnergy")),
      mantissaPrecision_(iConfig.getParameter<int>("mantissaPrecision")) {
  produces<Run3ScoutingCaloRecHitCollection>();
}

template<typename T>
void HLTScoutingCaloRecHitProducerT<T>::produce(edm::StreamID sid, edm::Event& iEvent, edm::EventSetup const& setup) const {

  auto const& recoCaloRecHits = iEvent.get(recoCaloRecHitsToken_);

  auto run3ScoutCaloRecHits = std::make_unique<Run3ScoutingCaloRecHitCollection>();
  run3ScoutCaloRecHits->reserve(recoCaloRecHits.size());

  for (auto const& recoCaloRecHit : recoCaloRecHits) {

    if (recoCaloRecHit.energy() < minEnergy_) {
      continue;
    }

    run3ScoutCaloRecHits->emplace_back(
      recoCaloRecHit.detid().rawId(),
      MiniFloatConverter::reduceMantissaToNbitsRounding(recoCaloRecHit.energy(), mantissaPrecision_)
    );
  }

  iEvent.put(std::move(run3ScoutCaloRecHits));
}

template<typename T>
void HLTScoutingCaloRecHitProducerT<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltCaloRecHits"));
  desc.add<double>("minEnergy", -1)->setComment("Minimum energy of the calorimeter RecHit in GeV");
  desc.add<int>("mantissaPrecision", 10)->setComment("default of 10 corresponds to float16 (change to 23 for float32)");
  descriptions.addWithDefaultLabel(desc);
}
