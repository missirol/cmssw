#ifndef CommonTools_Puppi_PuppiProducer_h
#define CommonTools_Puppi_PuppiProducer_h

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/PileupAlgos/interface/PuppiContainer.h"

class PuppiProducer : public edm::stream::EDProducer<> {
public:
  explicit PuppiProducer(edm::ParameterSet const&);
  ~PuppiProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  typedef math::XYZTLorentzVector LorentzVector;
  typedef std::vector<LorentzVector> LorentzVectorCollection;

private:
  void produce(edm::Event&, edm::EventSetup const&) override;

  edm::EDGetTokenT<reco::CandidateView> tokenPFCandidates_;
  edm::EDGetTokenT<reco::VertexCollection> tokenVertices_;
  edm::EDGetTokenT<reco::PFCandidateCollection> tokenPuppiCandidates_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> tokenPackedPuppiCandidates_;
  edm::EDPutTokenT<edm::ValueMap<float>> ptokenPupOut_;
  edm::EDPutTokenT<edm::ValueMap<LorentzVector>> ptokenP4PupOut_;
  edm::EDPutTokenT<edm::ValueMap<reco::CandidatePtr>> ptokenValues_;
  edm::EDPutTokenT<pat::PackedCandidateCollection> ptokenPackedPuppiCandidates_;
  edm::EDPutTokenT<reco::PFCandidateCollection> ptokenPuppiCandidates_;
  edm::EDPutTokenT<int> ptokenNalgos_;
  edm::EDPutTokenT<std::vector<double>> ptokenRawAlphas_;
  edm::EDPutTokenT<std::vector<double>> ptokenAlphas_;
  edm::EDPutTokenT<std::vector<double>> ptokenAlphasMed_;
  edm::EDPutTokenT<std::vector<double>> ptokenAlphasRms_;
  bool fPuppiDiagnostics;
  bool fPuppiNoLep;
  bool fUseFromPVLooseTight;
  bool fUseDZ;
  double fDZCut;
  double fEtaMinUseDZ;
  double fPtMaxCharged;
  double fEtaMaxCharged;
  double fPtMaxPhotons;
  double fEtaMaxPhotons;
  bool fUseExistingWeights;
  bool fClonePackedCands;
  int fVtxNdofCut;
  double fVtxZCut;
  std::unique_ptr<PuppiContainer> fPuppiContainer;
};

#endif
