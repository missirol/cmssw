#include <algorithm>
#include <memory>
#include <utility>

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "L1TriggerScouting/Utilities/interface/convertL1TObjectToL1ScoutingObjectT.h"

template <class T1, class T2>
class L1ScoutingObjectBXVecConverterT : public edm::global::EDProducer<> {
public:
  explicit L1ScoutingObjectBXVecConverterT(edm::ParameterSet const&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  edm::EDGetTokenT<BXVector<T1>> const srcToken_;
  int const bxMin_;
  int const bxMax_;
};

template <class T1, class T2>
L1ScoutingObjectBXVecConverterT<T1, T2>::L1ScoutingObjectBXVecConverterT(edm::ParameterSet const& iConfig)
    : srcToken_{consumes(iConfig.getParameter<edm::InputTag>("src"))},
      bxMin_{iConfig.getParameter<int>("bxMin")},
      bxMax_{iConfig.getParameter<int>("bxMax")} {
  produces<BXVector<T2>>();
}

template <class T1, class T2>
void L1ScoutingObjectBXVecConverterT<T1, T2>::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  auto const& input{iEvent.get(srcToken_)};

  auto const bxMin{std::max(bxMin_, input.getFirstBX())};
  auto const bxMax{std::min(bxMax_, input.getLastBX())};

  auto output{std::make_unique<BXVector<T2>>(0, bxMin, bxMax)};

  for (auto bx = bxMin; bx <= bxMax; ++bx) {
    auto const nInput{input.size(bx)};
    for (auto idx = 0u; idx < nInput; ++idx) {
      output->push_back(bx, convertL1TObjectToL1ScoutingObjectT<T1, T2>(input.at(bx, idx)));
    }
  }

  iEvent.put(std::move(output));
}

template <class T1, class T2>
void L1ScoutingObjectBXVecConverterT<T1, T2>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src")->setComment("Input collection (BXVector of L1T objects)");
  desc.add<int>("bxMin", -99)->setComment("Min BX (inclusive)");
  desc.add<int>("bxMax", 99)->setComment("Max BX (inclusive)");

  descriptions.addWithDefaultLabel(desc);
}
