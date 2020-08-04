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
  for (auto const& algo_i : lAlgos) {
    fPuppiAlgo.push_back(PuppiAlgo(algo_i));
  }
}

PuppiContainer::~PuppiContainer() {}

void PuppiContainer::initialize(std::vector<PuppiCandidate> const& iPuppiCands) {
  //Clear everything
  fWeights.clear();
  fVals.clear();
  fRawAlphas.clear();
  fAlphaMed.clear();
  fAlphaRMS.clear();
  fNPV = 1.;
  fCands.clear();
  fCands.reserve(iPuppiCands.size());
  fCandsForVar.clear();
  fCandsForVar.reserve(iPuppiCands.size());
  fCandsForVarChargedPV.clear();
  fCandsForVarChargedPV.reserve(iPuppiCands.size());
  for (auto const& pCand : iPuppiCands) {
    fCands.emplace_back(pCand);
    //charged candidates associated to PV
    if (pCand.id == 3)
      continue;
    fCandsForVar.emplace_back(pCand);
    if (pCand.id == 1)
      fCandsForVarChargedPV.emplace_back(pCand);
  }
}

float PuppiContainer::goodVar(PuppiCandidate const& iPuppiCand_0,
                              std::vector<PuppiCandidate> const& iPuppiCands,
                              int const iId,
                              float const iRCone) const {
  if (iId == -1)
    return 1.;

  float const r2(iRCone * iRCone);
  float var((iId == 1) ? iPuppiCand_0.pt : 0.);

  for (auto const& cand : iPuppiCands) {
    if (std::abs(cand.eta - iPuppiCand_0.eta) < iRCone) {
      auto const dr2(reco::deltaR2(cand.eta, cand.phi, iPuppiCand_0.eta, iPuppiCand_0.phi));
      if ((dr2 < r2) and (dr2 > 0.0001)) {
        auto const pt(cand.pt);
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
                               std::vector<PuppiCandidate> const& iPuppiCands,
                               std::vector<PuppiCandidate> const& iPuppiCandsForVar,
                               std::vector<PuppiCandidate> const& iPuppiCandsForVarChargedPV) {
  for (uint i0 = 0; i0 < iPuppiCands.size(); i0++) {
    auto const& iCand(iPuppiCands.at(i0));
    //Calculate the Puppi Algo to use
    auto const pPupId = getPuppiId(iCand.pt, iCand.eta);
    if (pPupId == -1 || fPuppiAlgo[pPupId].numAlgos() <= iOpt) {
      fVals.push_back(-1);
      continue;
    }

    //Get the Puppi Sub Algo (given iteration)
    auto const pAlgo = fPuppiAlgo[pPupId].algoId(iOpt);
    auto const pCharged = fPuppiAlgo[pPupId].isCharged(iOpt);
    auto const pCone = fPuppiAlgo[pPupId].coneSize(iOpt);
    //Compute the Puppi Metric
    float pVal(-1.);
    // calculate goodVar only for candidates that (1) will not be assigned a predefined weight (e.g 0, 1),
    // or (2) are required for computations inside puppi-algos (see call to PuppiAlgo::add below)
    if (((not(fApplyCHS and ((iCand.id == 1) or (iCand.id == 2)))) and (iCand.id != 3)) or
        ((std::abs(iCand.eta) < fPuppiAlgo[pPupId].etaMaxExtrap()) and ((iCand.id == 1) or (iCand.id == 2)))) {
      if (pCharged)
        pVal = goodVar(iCand, iPuppiCandsForVarChargedPV, pAlgo, pCone);
      else
        pVal = goodVar(iCand, iPuppiCandsForVar, pAlgo, pCone);
    }
    fVals.push_back(pVal);

    if (!edm::isFinite(pVal)) {
      LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iCand.pt << " -- " << iCand.eta;
      continue;
    }

    // // fPuppiAlgo[pPupId].add(iCand,pVal,iOpt);
    //code added by Nhan, now instead for every algorithm give it all the candidates
    for (int i1 = 0; i1 < fNAlgos; i1++) {
      // skip cands outside of algo's etaMaxExtrap,
      // as they would anyway be ignored inside PuppiAlgo::add
      if (not((std::abs(iCand.eta) < fPuppiAlgo[i1].etaMaxExtrap()) and ((iCand.id == 1) or (iCand.id == 2))))
        continue;

      auto curVal(pVal);
      if (i1 != pPupId) { //else, no need to repeat the computation
        auto const pAlgo = fPuppiAlgo[i1].algoId(iOpt);
        auto const pCharged = fPuppiAlgo[i1].isCharged(iOpt);
        auto const pCone = fPuppiAlgo[i1].coneSize(iOpt);
        if (pCharged)
          curVal = goodVar(iCand, iPuppiCandsForVarChargedPV, pAlgo, pCone);
        else
          curVal = goodVar(iCand, iPuppiCandsForVar, pAlgo, pCone);
      }

      fPuppiAlgo[i1].add(iCand, curVal, iOpt);
    }
  }

  for (int i0 = 0; i0 < fNAlgos; i0++)
    fPuppiAlgo[i0].computeMedRMS(iOpt);
}

void PuppiContainer::getRawAlphas(int const iOpt,
                                  std::vector<PuppiCandidate> const& iPuppiCands,
                                  std::vector<PuppiCandidate> const& iPuppiCandsForVar,
                                  std::vector<PuppiCandidate> const& iPuppiCandsForVarChargedPV) {
  for (int j0 = 0; j0 < fNAlgos; j0++) {
    for (auto const& iCand : iPuppiCands) {
      //Get the Puppi Sub Algo (given iteration)
      auto const pAlgo = fPuppiAlgo[j0].algoId(iOpt);
      auto const pCharged = fPuppiAlgo[j0].isCharged(iOpt);
      auto const pCone = fPuppiAlgo[j0].coneSize(iOpt);
      //Compute the Puppi Metric
      float pVal = -1.f;
      if (pCharged)
        pVal = goodVar(iCand, iPuppiCandsForVarChargedPV, pAlgo, pCone);
      else
        pVal = goodVar(iCand, iPuppiCandsForVar, pAlgo, pCone);
      fRawAlphas.push_back(pVal);
      if (!edm::isFinite(pVal)) {
        edm::LogWarning("NotFound") << "====> Value is Nan " << pVal << " == " << iCand.pt << " -- " << iCand.eta;
      }
    }
  }
}

int PuppiContainer::getPuppiId(float const iPt, float const iEta) {
  int lId = -1;
  for (int i0 = 0; i0 < fNAlgos; i0++) {
    auto const nEtaBinsPerAlgo = fPuppiAlgo[i0].etaBins();
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
  float lProbPU = 1.f - lProbLV;
  if (lProbPU <= 0)
    lProbPU = 1e-16;  //Quick Trick to through out infs
  if (lProbPU >= 0)
    lProbPU = 1 - 1e-16;  //Ditto
  float lChi2PU = TMath::ChisquareQuantile(lProbPU, 1);
  lChi2PU *= lChi2PU;
  return lChi2PU;
}

std::vector<float> const& PuppiContainer::puppiWeights() {
  auto const lNCands = fCands.size();

  fWeights.clear();
  fAlphaMed.clear();
  fAlphaRMS.clear();

  fWeights.reserve(lNCands);
  fAlphaMed.reserve(lNCands);
  fAlphaRMS.reserve(lNCands);

  for (int i0 = 0; i0 < fNAlgos; i0++)
    fPuppiAlgo[i0].reset();

  int lNMaxAlgo = 1;
  for (int i0 = 0; i0 < fNAlgos; i0++)
    lNMaxAlgo = std::max(fPuppiAlgo[i0].numAlgos(), lNMaxAlgo);

  //Run through all compute mean and RMS
  fVals.clear();
  fVals.reserve(lNCands * lNMaxAlgo);
  for (int i0 = 0; i0 < lNMaxAlgo; i0++)
    getRMSAvg(i0, fCands, fCandsForVar, fCandsForVarChargedPV);

  if (fPuppiDiagnostics)
    getRawAlphas(0, fCands, fCandsForVar, fCandsForVarChargedPV);

  for (uint i0 = 0; i0 < lNCands; i0++) {
    //Get the Puppi Id and, if ill defined, move on
    auto const& pCand(fCands.at(i0));
    auto const pPupId = getPuppiId(pCand.pt, pCand.eta);
    if (pPupId == -1) {
      fWeights.push_back(0);
      fAlphaMed.push_back(-10);
      fAlphaRMS.push_back(-10);
      continue;
    }

    float pWeight(1.);
    //Apply weight of 1 for leptons if puppiNoLep
    if (pCand.id == 3)
      pWeight = 1;
    //Apply the CHS weights
    else if ((pCand.id == 1) && fApplyCHS)
      pWeight = 1;
    else if ((pCand.id == 2) && fApplyCHS)
      pWeight = 0;
    else {
      // fill the p-values
      float pChi2 = 0;
      if (fUseExp) {
        //Compute an Experimental Puppi Weight with delta Z info (very simple example)
        pChi2 = getChi2FromdZ(pCand.dZ);
        //Now make sure Neutrals are not set
        if ((std::abs(pCand.pdgId) == 22) || (std::abs(pCand.pdgId) == 130))
          pChi2 = 0;
      }
      std::vector<float> pVals;
      auto const lNAlgos = fPuppiAlgo[pPupId].numAlgos();
      pVals.reserve(lNAlgos);
      for (int i1 = 0; i1 < lNAlgos; i1++) {
        pVals.push_back(fVals[lNCands * i1 + i0]);
      }
      //Fill and compute the PuppiWeight
      pWeight = fPuppiAlgo[pPupId].compute(pVals, pChi2);
    }

    //Basic Weight Checks
    if (not edm::isFinite(pWeight)) {
      pWeight = 0.0;
      LogDebug("PuppiWeightError") << "====> Weight is nan : " << pWeight << " : pt " << pCand.pt
                                   << " -- eta : " << pCand.eta << " -- Value" << fVals[i0] << " -- id :  " << pCand.id
                                   << " --  NAlgos: " << fPuppiAlgo[pPupId].numAlgos();
    }

    //Basic Cuts
    if (((pWeight * pCand.pt) < fPuppiAlgo[pPupId].neutralPt(fNPV)) && (pCand.id == 0))
      pWeight = 0;  //threshold cut on the neutral Pt

    // Protect high pT photons (important for gamma to hadronic recoil balance)
    if ((fPtMaxPhotons > 0) && (pCand.pdgId == 22) && (std::abs(pCand.eta) < fEtaMaxPhotons) &&
        (pCand.pt > fPtMaxPhotons))
      pWeight = 1.;
    // Protect high pT neutrals
    else if ((fPtMaxNeutrals > 0) && (pCand.id == 0))
      pWeight =
          std::clamp((pCand.pt - fPtMaxNeutralsStartSlope) / (fPtMaxNeutrals - fPtMaxNeutralsStartSlope), pWeight, 1.f);

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
