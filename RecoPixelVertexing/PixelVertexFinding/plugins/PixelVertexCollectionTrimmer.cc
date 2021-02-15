#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoPixelVertexing/PixelVertexFinding/interface/PVClusterComparer.h"

class PixelVertexCollectionTrimmer : public edm::stream::EDProducer<> {
public:
  explicit PixelVertexCollectionTrimmer(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<reco::VertexCollection> const vtxToken_;
  uint const maxVtx_;
  double const fractionSumPt2_;
  double const minSumPt2_;

  std::unique_ptr<PVClusterComparer> pvComparer_;
};

PixelVertexCollectionTrimmer::PixelVertexCollectionTrimmer(const edm::ParameterSet& iConfig)
    : vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      maxVtx_(iConfig.getParameter<uint>("maxVtx")),
      fractionSumPt2_(iConfig.getParameter<double>("fractionSumPt2")),
      minSumPt2_(iConfig.getParameter<double>("minSumPt2")) {
  auto const PVcomparerPSet = iConfig.getParameter<edm::ParameterSet>("PVcomparer");
  auto const track_pt_min = PVcomparerPSet.getParameter<double>("track_pt_min");
  auto const track_pt_max = PVcomparerPSet.getParameter<double>("track_pt_max");
  auto const track_chi2_max = PVcomparerPSet.getParameter<double>("track_chi2_max");
  auto const track_prob_min = PVcomparerPSet.getParameter<double>("track_prob_min");

  pvComparer_.reset(new PVClusterComparer(track_pt_min, track_pt_max, track_chi2_max, track_prob_min));

  produces<reco::VertexCollection>();
}

void PixelVertexCollectionTrimmer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto vtxs_trim = std::make_unique<reco::VertexCollection>();

  auto const vtxs = iEvent.getHandle(vtxToken_);

  if (vtxs->empty()) {
    edm::LogWarning("Input") << "Input collection of vertices is empty. Output collection will be empty.";
  } else {
    std::vector<double> foms(vtxs->size());
    for (size_t idx = 0; idx < vtxs->size(); ++idx)
      foms.at(idx) = pvComparer_->pTSquaredSum(vtxs->at(idx));

    std::vector<size_t> sorted_idxs(vtxs->size());
    std::iota(sorted_idxs.begin(), sorted_idxs.end(), 0);
    std::sort(
        sorted_idxs.begin(), sorted_idxs.end(), [&](size_t const i1, size_t const i2) { return foms[i1] > foms[i2]; });

    auto const min_fom = std::max(foms.at(sorted_idxs.at(0)) * fractionSumPt2_, minSumPt2_);

    vtxs_trim->reserve(vtxs->size());
    for (auto const idx : sorted_idxs) {
      if (vtxs_trim->size() >= maxVtx_)
        break;
      if (foms.at(idx) >= min_fom)
        vtxs_trim->emplace_back(vtxs->at(idx));
    }
  }

  iEvent.put(std::move(vtxs_trim));
}

void PixelVertexCollectionTrimmer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag(""))->setComment("input (pixel) vertex collection");
  desc.add<unsigned int>("maxVtx", 100)->setComment("max output collection size (number of accepted vertices)");
  desc.add<double>("fractionSumPt2", 0.3)->setComment("threshold on sumPt2 fraction of the leading vertex");
  desc.add<double>("minSumPt2", 0.)->setComment("min sumPt2");
  edm::ParameterSetDescription PVcomparerPSet;
  PVcomparerPSet.add<double>("track_pt_min", 1.0)->setComment("min track p_T");
  PVcomparerPSet.add<double>("track_pt_max", 10.0)->setComment("max track p_T");
  PVcomparerPSet.add<double>("track_chi2_max", 99999.)->setComment("max track chi2");
  PVcomparerPSet.add<double>("track_prob_min", -1.)->setComment("min track prob");
  desc.add<edm::ParameterSetDescription>("PVcomparer", PVcomparerPSet)
      ->setComment("from RecoPixelVertexing/PixelVertexFinding/python/PVClusterComparer_cfi.py");
  descriptions.add("hltPixelVertexCollectionTrimmer", desc);
}

DEFINE_FWK_MODULE(PixelVertexCollectionTrimmer);
