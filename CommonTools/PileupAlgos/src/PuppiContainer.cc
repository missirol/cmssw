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
  fWeights.clear();
  fVals.clear();
  fRawAlphas.clear();
  fAlphaMed.clear();
  fAlphaRMS.clear();
  fNPV = 1.;

  fCands = iPuppiCandidates;
  fCandIndicesSum.clear();
  fCandIndicesSum.reserve(iPuppiCandidates.size());
  fCandIndicesSumChargedPV.clear();
  fCandIndicesSumChargedPV.reserve(iPuppiCandidates.size());
  for (uint idx=0; idx<fCands.size(); ++idx) {
    if (fCands.at(idx).id == 3)
      continue;
    fCandIndicesSum.emplace_back(idx);
    // indices of charged candidates associated to LV
    if (fCands.at(idx).id == 1)
      fCandIndicesSumChargedPV.emplace_back(idx);
  }
}

double PuppiContainer::goodVar(std::vector<PuppiCandidate> const& iCands,
			       uint const iCandIdx,
			       std::vector<uint> const& iCandIndices,
                               int const iId,
                               double const iRCone) const {
  if (iId == -1)
    return 1.;

  auto const& cand_0(iCands.at(iCandIdx));

  double const r2(iRCone * iRCone);
  double var((iId == 1) ? cand_0.pt : 0.);

  for (auto const candIdx_j : iCandIndices) {
    if(candIdx_j == iCandIdx)
      continue;
    auto const cand(iCands.at(candIdx_j));
    if (std::abs(cand.eta - cand_0.eta) < iRCone) {
      auto const dr2(reco::deltaR2(cand.eta, cand.phi, cand_0.eta, cand_0.phi));
      if ((dr2 < r2) and (dr2 > 0.0001)) { //!!
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
                               std::vector<PuppiCandidate> const& iCands,
                               std::vector<uint> const& iCandIndicesSum,
                               std::vector<uint> const& iCandIndicesSumChargedPV) {
  for (uint i0 = 0; i0 < iCands.size(); i0++) {
    auto const& iCand(iCands.at(i0));
    //Calculate the Puppi Algo to use
    int pPupId = getPuppiId(iCand.pt, iCand.eta);
    if (pPupId == -1 || fPuppiAlgo[pPupId].numAlgos() <= iOpt) {
      fVals.push_back(-1);
      continue;
    }

    //Get the Puppi Sub Algo (given iteration)
    int pAlgo = fPuppiAlgo[pPupId].algoId(iOpt);
    bool pCharged = fPuppiAlgo[pPupId].isCharged(iOpt);
    double pCone = fPuppiAlgo[pPupId].coneSize(iOpt);
    //Compute the Puppi Metric
    double pVal(-1.);
    // calculate goodVar only for candidates that (1) will not be assigned a predefined weight (e.g. 0, 1),
    // or (2) are required for computations inside puppi-algos (see call to PuppiAlgo::add below)
    if (((not(fApplyCHS and ((iCand.id == 1) or (iCand.id == 2)))) and (iCand.id != 3)) or
        ((std::abs(iCand.eta) < fPuppiAlgo[pPupId].etaMaxExtrap()) and ((iCand.id == 1) or (iCand.id == 2)))) {
      if (pCharged)
        pVal = goodVar(iCands, i0, iCandIndicesSumChargedPV, pAlgo, pCone);
      else
        pVal = goodVar(iCands, i0, iCandIndicesSum, pAlgo, pCone);
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
      pAlgo = fPuppiAlgo[i1].algoId(iOpt);
      pCharged = fPuppiAlgo[i1].isCharged(iOpt);
      pCone = fPuppiAlgo[i1].coneSize(iOpt);
      double curVal = -1;
      if (i1 != pPupId) {
        if (pCharged)
          curVal = goodVar(iCands, i0, iCandIndicesSumChargedPV, pAlgo, pCone);
        else
          curVal = goodVar(iCands, i0, iCandIndicesSum, pAlgo, pCone);
      } else {  //no need to repeat the computation
        curVal = pVal;
      }

      fPuppiAlgo[i1].add(iCand, curVal, iOpt);
    }
  }

  for (int i0 = 0; i0 < fNAlgos; i0++)
    fPuppiAlgo[i0].computeMedRMS(iOpt);
}

void PuppiContainer::getRawAlphas(int const iOpt,
                                  std::vector<PuppiCandidate> const& iCands,
                                  std::vector<uint> const& iCandIndicesSum,
                                  std::vector<uint> const& iCandIndicesSumChargedPV) {
  for (int j0 = 0; j0 < fNAlgos; j0++) {
    for (uint candIdx=0; candIdx<iCands.size(); ++candIdx) {
      double pVal = -1;
      //Get the Puppi Sub Algo (given iteration)
      int pAlgo = fPuppiAlgo[j0].algoId(iOpt);
      bool pCharged = fPuppiAlgo[j0].isCharged(iOpt);
      double pCone = fPuppiAlgo[j0].coneSize(iOpt);
      //Compute the Puppi Metric
      if (pCharged)
        pVal = goodVar(iCands, candIdx, iCandIndicesSumChargedPV, pAlgo, pCone);
      else
        pVal = goodVar(iCands, candIdx, iCandIndicesSum, pAlgo, pCone);
      fRawAlphas.push_back(pVal);
      if (!edm::isFinite(pVal)) {
        LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iCands.at(candIdx).pt << " -- " << iCands.at(candIdx).eta;
        continue;
      }
    }
  }
}

int PuppiContainer::getPuppiId(double const iPt, double const iEta) {
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

double PuppiContainer::getChi2FromdZ(double const iDZ) const {
  //We need to obtain prob of PU + (1-Prob of LV)
  // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm  (its really more like 1mm)
  //double lProbLV = ROOT::Math::normal_cdf_c(std::abs(iDZ),0.2)*2.; //*2 is to do it double sided
  //Take iDZ to be corrected by sigma already
  double lProbLV = ROOT::Math::normal_cdf_c(std::abs(iDZ), 1.) * 2.;  //*2 is to do it double sided
  double lProbPU = 1 - lProbLV;
  if (lProbPU <= 0)
    lProbPU = 1e-16;  //Quick Trick to through out infs
  if (lProbPU >= 0)
    lProbPU = 1 - 1e-16;  //Ditto
  double lChi2PU = TMath::ChisquareQuantile(lProbPU, 1);
  lChi2PU *= lChi2PU;
  return lChi2PU;
}

std::vector<double> const& PuppiContainer::puppiWeights() {
  int const lNCands(fCands.size());

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
    getRMSAvg(i0, fCands, fCandIndicesSum, fCandIndicesSumChargedPV);

  if (fPuppiDiagnostics)
    getRawAlphas(0, fCands, fCandIndicesSum, fCandIndicesSumChargedPV);

  for (int i0 = 0; i0 < lNCands; i0++) {
    //Get the Puppi Id and, if ill defined, move on
    auto const& pCand(fCands.at(i0));
    int pPupId = getPuppiId(pCand.pt, pCand.eta);
    if (pPupId == -1) {
      fWeights.push_back(0);
      fAlphaMed.push_back(-10);
      fAlphaRMS.push_back(-10);
      continue;
    }

    double pWeight(1.);
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
      double pChi2 = 0;
      if (fUseExp) {
        //Compute an Experimental Puppi Weight with delta Z info (very simple example)
        pChi2 = getChi2FromdZ(pCand.dZ);
        //Now make sure Neutrals are not set
        if ((std::abs(pCand.pdgId) == 22) || (std::abs(pCand.pdgId) == 130))
          pChi2 = 0;
      }
      std::vector<double> pVals;
      int lNAlgos = fPuppiAlgo[pPupId].numAlgos();
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
          std::clamp((pCand.pt - fPtMaxNeutralsStartSlope) / (fPtMaxNeutrals - fPtMaxNeutralsStartSlope), pWeight, 1.);

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
