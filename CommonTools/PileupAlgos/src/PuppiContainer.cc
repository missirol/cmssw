#include <cmath>

#include "CommonTools/PileupAlgos/interface/PuppiContainer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/ProbFuncMathCore.h"
#include "TMath.h"

PuppiContainer::PuppiContainer(edm::ParameterSet const& iConfig) {
  fPuppiDiagnostics = iConfig.getParameter<bool>("puppiDiagnostics");
  fApplyCHS = iConfig.getParameter<bool>("applyCHS");
  fInvert = iConfig.getParameter<bool>("invertPuppi");
  fUseExp = iConfig.getParameter<bool>("useExp");
  fPuppiWeightCut = iConfig.getParameter<double>("MinPuppiWeight");
  fPtMaxPhotons = iConfig.getParameter<double>("PtMaxPhotons");
  fEtaMaxPhotons = iConfig.getParameter<double>("EtaMaxPhotons");
  fPtMaxNeutrals = iConfig.getParameter<double>("PtMaxNeutrals");
  fPtMaxNeutralsStartSlope = iConfig.getParameter<double>("PtMaxNeutralsStartSlope");
  std::vector<edm::ParameterSet> lAlgos = iConfig.getParameter<std::vector<edm::ParameterSet>>("algos");
  fNAlgos = lAlgos.size();
  fPuppiAlgo.reserve(lAlgos.size());
  for (uint i0 = 0; i0 < lAlgos.size(); i0++) {
    PuppiAlgo pPuppiConfig(lAlgos[i0]);
    fPuppiAlgo.push_back(pPuppiConfig);
  }
}

PuppiContainer::~PuppiContainer() {}

void PuppiContainer::initialize(std::vector<PuppiCandidate> const& iPuppiCandidates) {
  //Clear everything
  fPFParticles.clear();
  fChargedPV.clear();
  fWeights.clear();
  fVals.clear();
  fRawAlphas.clear();
  fAlphaMed.clear();
  fAlphaRMS.clear();
  fNPV = 1.;
  fPFParticles.reserve(iPuppiCandidates.size());
  fChargedPV.reserve(iPuppiCandidates.size());
  for (auto const& pParticle : iPuppiCandidates) {
    fPFParticles.emplace_back(pParticle);
    //charged particles associated to PV
    if (pParticle.id == 1)
      fChargedPV.emplace_back(pParticle);
  }
}

float PuppiContainer::goodVar(PuppiCandidate const& iPart0,
                               std::vector<PuppiCandidate> const& iParticles,
                               int const iId,
                               float const iRCone) const {
  if (iId == -1)
    return 1.;

  float const r2(iRCone * iRCone);
  float var((iId == 1) ? iPart0.pt : 0.);

  for (auto const& part : iParticles) {
    if (part.id == 3)
      continue;
    if (std::abs(part.eta - iPart0.eta) < iRCone) {
      auto const dr2(reco::deltaR2(part.eta, part.phi, iPart0.eta, iPart0.phi));
      if ((dr2 < r2) and (dr2 > 0.0001)) {
        auto const pt(part.pt);
        if (iId == 5)
          var += (pt * pt / dr2);
        else if (iId == 4)
          var += pt;
        else if (iId == 3)
          var += (1. / dr2);
        else if (iId == 2)
          var += (1. / dr2);
        else if (iId == 1)
          var += pt;
        else if (iId == 0)
          var += (pt / dr2);
      }
    }
  }

  if ((var != 0) and ((iId == 0) or (iId == 3) or (iId == 5)))
    var = log(var);

  return var;
}

//In fact takes the median not the average
void PuppiContainer::getRMSAvg(int const iOpt,
                               std::vector<PuppiCandidate> const& iParticles,
                               std::vector<PuppiCandidate> const& iChargedParticles) {
  for (uint i0 = 0; i0 < iParticles.size(); i0++) {
    auto const& iPart(iParticles.at(i0));
    //Calculate the Puppi Algo to use
    int pPupId = getPuppiId(iPart.pt, iPart.eta);
    if (pPupId == -1 || fPuppiAlgo[pPupId].numAlgos() <= iOpt) {
      fVals.push_back(-1);
      continue;
    }

    //Get the Puppi Sub Algo (given iteration)
    int pAlgo = fPuppiAlgo[pPupId].algoId(iOpt);
    bool pCharged = fPuppiAlgo[pPupId].isCharged(iOpt);
    float pCone = fPuppiAlgo[pPupId].coneSize(iOpt);
    //Compute the Puppi Metric
    float pVal(-1.);
    // calculate goodVar only for candidates that (1) will not be assigned a predefined weight (e.g 0, 1),
    // or (2) are required for computations inside puppi-algos (see call to PuppiAlgo::add below)
    if (((not(fApplyCHS and ((iPart.id == 1) or (iPart.id == 2)))) and (iPart.id != 3)) or
        ((std::abs(iPart.eta) < fPuppiAlgo[pPupId].etaMaxExtrap()) and ((iPart.id == 1) or (iPart.id == 2)))) {
      if (!pCharged)
        pVal = goodVar(iPart, iParticles, pAlgo, pCone);
      if (pCharged)
        pVal = goodVar(iPart, iChargedParticles, pAlgo, pCone);
    }
    fVals.push_back(pVal);

    if (!edm::isFinite(pVal)) {
      LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iPart.pt << " -- " << iPart.eta;
      continue;
    }

    // // fPuppiAlgo[pPupId].add(iPart,pVal,iOpt);
    //code added by Nhan, now instead for every algorithm give it all the particles
    for (int i1 = 0; i1 < fNAlgos; i1++) {
      // skip cands outside of algo's etaMaxExtrap,
      // as they would anyway be ignored inside PuppiAlgo::add
      if (not((std::abs(iPart.eta) < fPuppiAlgo[i1].etaMaxExtrap()) and ((iPart.id == 1) or (iPart.id == 2))))
        continue;

      pAlgo = fPuppiAlgo[i1].algoId(iOpt);
      pCharged = fPuppiAlgo[i1].isCharged(iOpt);
      pCone = fPuppiAlgo[i1].coneSize(iOpt);
      float curVal = -1;
      if (i1 != pPupId) {
        if (!pCharged)
          curVal = goodVar(iPart, iParticles, pAlgo, pCone);
        if (pCharged)
          curVal = goodVar(iPart, iChargedParticles, pAlgo, pCone);
      } else {  //no need to repeat the computation
        curVal = pVal;
      }

      fPuppiAlgo[i1].add(iPart, curVal, iOpt);
    }
  }

  for (int i0 = 0; i0 < fNAlgos; i0++)
    fPuppiAlgo[i0].computeMedRMS(iOpt);
}

void PuppiContainer::getRawAlphas(int const iOpt,
                                  std::vector<PuppiCandidate> const& iParticles,
                                  std::vector<PuppiCandidate> const& iChargedParticles) {
  for (int j0 = 0; j0 < fNAlgos; j0++) {
    for (auto const& iPart : iParticles) {
      float pVal = -1;
      //Get the Puppi Sub Algo (given iteration)
      int pAlgo = fPuppiAlgo[j0].algoId(iOpt);
      bool pCharged = fPuppiAlgo[j0].isCharged(iOpt);
      float pCone = fPuppiAlgo[j0].coneSize(iOpt);
      //Compute the Puppi Metric
      if (!pCharged)
        pVal = goodVar(iPart, iParticles, pAlgo, pCone);
      if (pCharged)
        pVal = goodVar(iPart, iChargedParticles, pAlgo, pCone);
      fRawAlphas.push_back(pVal);
      if (!edm::isFinite(pVal)) {
        LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iPart.pt << " -- " << iPart.eta;
        continue;
      }
    }
  }
}

int PuppiContainer::getPuppiId(float const iPt, float const iEta) {
  int lId = -1;
  for (int i0 = 0; i0 < fNAlgos; i0++) {
    int nEtaBinsPerAlgo = fPuppiAlgo[i0].etaBins();
    for (int i1 = 0; i1 < nEtaBinsPerAlgo; i1++) {
      if ((std::abs(iEta) > fPuppiAlgo[i0].etaMin(i1)) && (std::abs(iEta) < fPuppiAlgo[i0].etaMax(i1))) {
        fPuppiAlgo[i0].fixAlgoEtaBin(i1);
        if (iPt > fPuppiAlgo[i0].ptMin()) {
          lId = i0;
          break;
        }
      }
    }
  }

  return lId;
}

float PuppiContainer::getChi2FromdZ(float const iDZ) const {
  //We need to obtain prob of PU + (1-Prob of LV)
  // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm  (its really more like 1mm)
  //float lProbLV = ROOT::Math::normal_cdf_c(std::abs(iDZ),0.2)*2.; //*2 is to do it double sided
  //Take iDZ to be corrected by sigma already
  float lProbLV = ROOT::Math::normal_cdf_c(std::abs(iDZ), 1.) * 2.;  //*2 is to do it double sided
  float lProbPU = 1 - lProbLV;
  if (lProbPU <= 0)
    lProbPU = 1e-16;  //Quick Trick to through out infs
  if (lProbPU >= 0)
    lProbPU = 1 - 1e-16;  //Ditto
  float lChi2PU = TMath::ChisquareQuantile(lProbPU, 1);
  lChi2PU *= lChi2PU;
  return lChi2PU;
}

std::vector<float> const& PuppiContainer::puppiWeights() {
  int const lNParticles = fPFParticles.size();

  fWeights.clear();
  fAlphaMed.clear();
  fAlphaRMS.clear();

  fWeights.reserve(lNParticles);
  fAlphaMed.reserve(lNParticles);
  fAlphaRMS.reserve(lNParticles);

  for (int i0 = 0; i0 < fNAlgos; i0++)
    fPuppiAlgo[i0].reset();

  int lNMaxAlgo = 1;
  for (int i0 = 0; i0 < fNAlgos; i0++)
    lNMaxAlgo = std::max(fPuppiAlgo[i0].numAlgos(), lNMaxAlgo);

  //Run through all compute mean and RMS
  fVals.clear();
  fVals.reserve(lNParticles * lNMaxAlgo);
  for (int i0 = 0; i0 < lNMaxAlgo; i0++)
    getRMSAvg(i0, fPFParticles, fChargedPV);

  if (fPuppiDiagnostics)
    getRawAlphas(0, fPFParticles, fChargedPV);

  for (int i0 = 0; i0 < lNParticles; i0++) {
    //Get the Puppi Id and, if ill defined, move on
    auto const& rPart(fPFParticles.at(i0));
    int pPupId = getPuppiId(rPart.pt, rPart.eta);
    if (pPupId == -1) {
      fWeights.push_back(0);
      fAlphaMed.push_back(-10);
      fAlphaRMS.push_back(-10);
      continue;
    }

    float pWeight(1.);
    //Apply weight of 1 for leptons if puppiNoLep
    if (rPart.id == 3)
      pWeight = 1;
    //Apply the CHS weights
    else if ((rPart.id == 1) && fApplyCHS)
      pWeight = 1;
    else if ((rPart.id == 2) && fApplyCHS)
      pWeight = 0;
    else {
      // fill the p-values
      float pChi2 = 0;
      if (fUseExp) {
        //Compute an Experimental Puppi Weight with delta Z info (very simple example)
        pChi2 = getChi2FromdZ(rPart.dZ);
        //Now make sure Neutrals are not set
        if ((std::abs(rPart.pdgId) == 22) || (std::abs(rPart.pdgId) == 130))
          pChi2 = 0;
      }
      std::vector<float> pVals;
      int lNAlgos = fPuppiAlgo[pPupId].numAlgos();
      pVals.reserve(lNAlgos);
      for (int i1 = 0; i1 < lNAlgos; i1++) {
        pVals.push_back(fVals[lNParticles * i1 + i0]);
      }
      //Fill and compute the PuppiWeight
      pWeight = fPuppiAlgo[pPupId].compute(pVals, pChi2);
    }

    //Basic Weight Checks
    if (not edm::isFinite(pWeight)) {
      pWeight = 0.0;
      LogDebug("PuppiWeightError") << "====> Weight is nan : " << pWeight << " : pt " << rPart.pt
                                   << " -- eta : " << rPart.eta << " -- Value" << fVals[i0] << " -- id :  " << rPart.id
                                   << " --  NAlgos: " << fPuppiAlgo[pPupId].numAlgos();
    }

    //Basic Cuts
    if (((pWeight * rPart.pt) < fPuppiAlgo[pPupId].neutralPt(fNPV)) && (rPart.id == 0))
      pWeight = 0;  //threshold cut on the neutral Pt

    // Protect high pT photons (important for gamma to hadronic recoil balance)
    if ((fPtMaxPhotons > 0) && (rPart.pdgId == 22) && (std::abs(rPart.eta) < fEtaMaxPhotons) &&
        (rPart.pt > fPtMaxPhotons))
      pWeight = 1.;
    // Protect high pT neutrals
    else if ((fPtMaxNeutrals > 0) && (rPart.id == 0))
      pWeight =
          std::clamp((rPart.pt - fPtMaxNeutralsStartSlope) / (fPtMaxNeutrals - fPtMaxNeutralsStartSlope), pWeight, 1.f);

    // Eliminate the low Weight stuff
    if (pWeight < fPuppiWeightCut)
      pWeight = 0;

    if (fInvert)
      pWeight = 1. - pWeight;

    fWeights.push_back(pWeight);
    fAlphaMed.push_back(fPuppiAlgo[pPupId].median());
    fAlphaRMS.push_back(fPuppiAlgo[pPupId].rms());
  }

  return fWeights;
}
