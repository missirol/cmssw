#include "CommonTools/PileupAlgos/interface/PuppiContainer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/ProbFunc.h"
#include "TMath.h"
#include <iostream>
#include <cmath>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"

using namespace std;

PuppiContainer::PuppiContainer(const edm::ParameterSet &iConfig) {
  fPuppiDiagnostics = iConfig.getParameter<bool>("puppiDiagnostics");
  fApplyCHS = iConfig.getParameter<bool>("applyCHS");
  fInvert = iConfig.getParameter<bool>("invertPuppi");
  fUseExp = iConfig.getParameter<bool>("useExp");
  fPuppiWeightCut = iConfig.getParameter<double>("MinPuppiWeight");
  fPtMaxPhotons = iConfig.getParameter<double>("PtMaxPhotons");
  fEtaMaxPhotons = iConfig.getParameter<double>("EtaMaxPhotons");
  fPtMaxNeutrals = iConfig.getParameter<double>("PtMaxNeutrals");
  fPtMaxNeutralsStartSlope = iConfig.getParameter<double>("PtMaxNeutralsStartSlope");
  std::vector<edm::ParameterSet> lAlgos = iConfig.getParameter<std::vector<edm::ParameterSet> >("algos");
  fNAlgos = lAlgos.size();
  for (unsigned int i0 = 0; i0 < lAlgos.size(); i0++) {
    PuppiAlgo pPuppiConfig(lAlgos[i0]);
    fPuppiAlgo.push_back(pPuppiConfig);
  }
}

PuppiContainer::~PuppiContainer() {}

//template<unsigned int N>
//float PuppiContainer::goodVar(PuppiCandidates<N> const &candidates,
//                              uint const candIdx,
//                              bool const useLVCands,
//                              int const iId,
//                              float const R) {
//  if (iId == -1)
//    return 1.f;
//
//  float const r2 = R * R;
//  float var = 0.f;
//
//  for (uint idx=0; idx<candidates.size(); ++idx) {
//    if(idx == candIdx or candidates.id[candIdx] == 3 or (useLVCands and candidates.id[candIdx] == 1))
//      continue;
//
//    if (std::abs(candidates.rapidity[idx] - candidates.rapidity[candIdx]) < R) {
//      auto const dr2y = reco::deltaR2(candidates.rapidity[idx], candidates.phi[idx], candidates.rapidity[candIdx], candidates.phi[candIdx]);
//      if (dr2y < r2) {
//        auto const dr2 = reco::deltaR2(candidates.eta[idx], candidates.phi[idx], candidates.eta[candIdx], candidates.phi[candIdx]);
//        auto const pt = candidates.pt[idx];
//        if (iId == 5)
//          var += (pt * pt / dr2);
//        else if (iId == 4)
//          var += pt;
//        else if (iId == 3)
//          var += (1. / dr2);
//        else if (iId == 2)
//          var += (1. / dr2);
//        else if (iId == 1)
//          var += pt;
//        else if (iId == 0)
//          var += (pt / dr2);
//      }
//    }
//  }
//
//  if ((var != 0.) and ((iId == 0) or (iId == 3) or (iId == 5)))
//    var = log(var);
//  else if (iId == 1)
//    var += candidates.pt[candIdx];
//
//  return var;
//}

//In fact takes the median not the average
//void PuppiContainer::getRMSAvg(int iOpt, PuppiCandidates const &iConstits) {
//  for (unsigned int i0 = 0; i0 < iConstits.size(); i0++) {
//    //Calculate the Puppi Algo to use
//    int pPupId = getPuppiId(iConstits[i0].pt, iConstits[i0].eta);
//    if (pPupId == -1 || fPuppiAlgo[pPupId].numAlgos() <= iOpt) {
//      fVals.push_back(-1);
//      continue;
//    }
//    //Get the Puppi Sub Algo (given iteration)
//    int pAlgo = fPuppiAlgo[pPupId].algoId(iOpt);
//    bool pCharged = fPuppiAlgo[pPupId].isCharged(iOpt);
//    float pCone = fPuppiAlgo[pPupId].coneSize(iOpt);
//    // compute the Puppi metric:
//    //  - calculate goodVar only for candidates that (1) will not be assigned a predefined weight (e.g 0, 1),
//    //    or (2) are required for computations inside puppi-algos (see call to PuppiAlgo::add below)
//    float pVal = -1;
//    bool const getsDefaultWgtIfApplyCHS = iConstits[i0].id == 1 or iConstits[i0].id == 2;
//    if (not((fApplyCHS and getsDefaultWgtIfApplyCHS) or iConstits[i0].id == 3) or
//        (std::abs(iConstits[i0].eta) < fPuppiAlgo[pPupId].etaMaxExtrap() and getsDefaultWgtIfApplyCHS)) {
//      pVal = goodVar(iConstits, i0, pCharged, pAlgo, pCone);
//    }
//    fVals.push_back(pVal);
//
//    if (!edm::isFinite(pVal)) {
//      LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iConstits[i0].pt << " -- " << iConstits[i0].eta
//                           << endl;
//      continue;
//    }
//
//    // code added by Nhan: now instead for every algorithm give it all the particles
//    for (int i1 = 0; i1 < fNAlgos; i1++) {
//      // skip cands outside of algo's etaMaxExtrap, as they would anyway be ignored inside PuppiAlgo::add (see end of the block)
//      if (not(std::abs(iConstits[i0].eta) < fPuppiAlgo[i1].etaMaxExtrap() and getsDefaultWgtIfApplyCHS))
//        continue;
//
//      auto curVal = pVal;
//      // recompute goodVar if algo has changed
//      if (i1 != pPupId) {
//        pAlgo = fPuppiAlgo[i1].algoId(iOpt);
//        pCharged = fPuppiAlgo[i1].isCharged(iOpt);
//        pCone = fPuppiAlgo[i1].coneSize(iOpt);
//        curVal = goodVar(iConstits, i0, pCharged, pAlgo, pCone);
//      }
//
//      fPuppiAlgo[i1].add(curVal, iOpt, iConstits.pt[i0], iConstits.eta[i0], iConstits.id[i0]);
//    }
//  }
//
//  for (int i0 = 0; i0 < fNAlgos; i0++)
//    fPuppiAlgo[i0].computeMedRMS(iOpt);
//}

//In fact takes the median not the average
//void PuppiContainer::getRawAlphas(int iOpt, PuppiCandidates const &iConstits) {
//  for (int j0 = 0; j0 < fNAlgos; j0++) {
//    for (unsigned int i0 = 0; i0 < iConstits.size(); i0++) {
//      //Get the Puppi Sub Algo (given iteration)
//      int pAlgo = fPuppiAlgo[j0].algoId(iOpt);
//      bool pCharged = fPuppiAlgo[j0].isCharged(iOpt);
//      float pCone = fPuppiAlgo[j0].coneSize(iOpt);
//      //Compute the Puppi Metric
//      float const pVal = goodVar(iConstits, i0, pCharged, pAlgo, pCone);
//      fRawAlphas.push_back(pVal);
//      if (!edm::isFinite(pVal)) {
//        LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iConstits[i0].pt << " -- "
//                             << iConstits[i0].eta << endl;
//        continue;
//      }
//    }
//  }
//}

int PuppiContainer::getPuppiId(float iPt, float iEta) {
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
  //if(lId == -1) std::cerr << "Error : Full fiducial range is not defined " << std::endl;
  return lId;
}
float PuppiContainer::getChi2FromdZ(float iDZ) {
  //We need to obtain prob of PU + (1-Prob of LV)
  // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm  (its really more like 1mm)
  //float lProbLV = ROOT::Math::normal_cdf_c(std::abs(iDZ),0.2)*2.; //*2 is to do it double sided
  //Take iDZ to be corrected by sigma already
  float lProbLV = ROOT::Math::normal_cdf_c(std::abs(iDZ), 1.) * 2.f;  //*2 is to do it double sided
  float lProbPU = 1 - lProbLV;
  if (lProbPU <= 0)
    lProbPU = 1e-16;  //Quick Trick to through out infs
  if (lProbPU >= 0)
    lProbPU = 1 - 1e-16;  //Ditto
  float lChi2PU = TMath::ChisquareQuantile(lProbPU, 1);
  lChi2PU *= lChi2PU;
  return lChi2PU;
}

int PuppiContainer::fillValues(std::vector<RecoObj> const& recoObjs, int const numberOfPVs) {
  //Clear everything
  fWeights.resize(0);
  fVals.resize(0);
  fRawAlphas.resize(0);
  fAlphaMed.resize(0);
  fAlphaRMS.resize(0);

  uint const lNParticles = recoObjs.size();

  PuppiCandidates<lNParticles> puppiCandidates;
  for (uint idx=0; idx<recoObjs.size(); ++idx) {
    auto const& rParticle = recoObjs.at(idx);
    puppiCandidates.id[idx] = rParticle.id;
    if (edm::isFinite(rParticle.rapidity)) {
      puppiCandidates.pt[idx] = rParticle.pt;
      puppiCandidates.eta[idx] = rParticle.eta;
      puppiCandidates.rapidity[idx] = rParticle.rapidity;
      puppiCandidates.phi[idx] = rParticle.phi;
    } else {
      puppiCandidates.pt[idx] = 0.f;
      puppiCandidates.eta[idx] = 99.f;
      puppiCandidates.rapidity[idx] = 99.f;
      puppiCandidates.phi[idx] = 0.f;
    }
  }

  fWeights.reserve(lNParticles);
  fVals.reserve(lNParticles);
  for (int i0 = 0; i0 < fNAlgos; i0++)
    fPuppiAlgo[i0].reset();

  int lNMaxAlgo = 1;
  for (int i0 = 0; i0 < fNAlgos; i0++)
    lNMaxAlgo = std::max(fPuppiAlgo[i0].numAlgos(), lNMaxAlgo);
  //Run through all compute mean and RMS
  for (int i0 = 0; i0 < lNMaxAlgo; i0++) {
    getRMSAvg<lNParticles>(i0, puppiCandidates);
  }
  if (fPuppiDiagnostics)
    getRawAlphas<lNParticles>(0, puppiCandidates);

  std::vector<float> pVals;
  pVals.reserve(lNParticles);
  for (uint i0 = 0; i0 < lNParticles; i0++) {
    auto const& rParticle = recoObjs.at(i0);
    //Refresh
    pVals.clear();
    float pWeight = 1;
    //Get the Puppi Id and if ill defined move on
    int pPupId = getPuppiId(puppiCandidates.pt[i0], puppiCandidates.eta[i0]);
    if (pPupId == -1) {
      fWeights.push_back(0);
      fAlphaMed.push_back(-10);
      fAlphaRMS.push_back(-10);
      continue;
    }

    // fill the p-values
    float pChi2 = 0;
    if (fUseExp) {
      //Compute an Experimental Puppi Weight with delta Z info (very simple example)
      pChi2 = getChi2FromdZ(rParticle.dZ);
      //Now make sure Neutrals are not set
      if ((std::abs(rParticle.pdgId) == 22) || (std::abs(rParticle.pdgId) == 130))
        pChi2 = 0;
    }
    //Fill and compute the PuppiWeight
    int lNAlgos = fPuppiAlgo[pPupId].numAlgos();
    for (int i1 = 0; i1 < lNAlgos; i1++)
      pVals.push_back(fVals[lNParticles * i1 + i0]);

    pWeight = fPuppiAlgo[pPupId].compute(pVals, pChi2);
    //Apply the CHS weights
    if (puppiCandidates.id[i0] == 1 && fApplyCHS)
      pWeight = 1;
    if (puppiCandidates.id[i0] == 2 && fApplyCHS)
      pWeight = 0;
    //Apply weight of 1 for leptons if puppiNoLep
    if (puppiCandidates.id[i0] == 3)
      pWeight = 1;
    //Basic Weight Checks
    if (!edm::isFinite(pWeight)) {
      pWeight = 0.0;
      LogDebug("PuppiWeightError") << "====> Weight is nan : " << pWeight << " : pt " << puppiCandidates.pt[i0]
                                   << " -- eta : " << puppiCandidates.eta[i0] << " -- Value" << fVals[i0]
                                   << " -- id :  " << puppiCandidates.id[i0] << " --  NAlgos: " << lNAlgos << std::endl;
    }
    //Basic Cuts
    if (pWeight * rParticle.pt < fPuppiAlgo[pPupId].neutralPt(numberOfPVs) && puppiCandidates.id[i0] == 0)
      pWeight = 0;  //threshold cut on the neutral Pt
    // Protect high pT photons (important for gamma to hadronic recoil balance)
    if (fPtMaxPhotons > 0 && rParticle.pdgId == 22 && std::abs(rParticle.eta) < fEtaMaxPhotons && rParticle.pt > fPtMaxPhotons)
      pWeight = 1.f;
    // Protect high pT neutrals
    else if ((fPtMaxNeutrals > 0) && (puppiCandidates.id[i0] == 0))
      pWeight = std::clamp((rParticle.pt-fPtMaxNeutralsStartSlope)/(fPtMaxNeutrals-fPtMaxNeutralsStartSlope), pWeight, 1.f);

    if (pWeight < fPuppiWeightCut)
      pWeight = 0;  //==> Elminate the low Weight stuff
    if (fInvert)
      pWeight = 1. - pWeight;

    fWeights.push_back(pWeight);
    fAlphaMed.push_back(fPuppiAlgo[pPupId].median());
    fAlphaRMS.push_back(fPuppiAlgo[pPupId].rms());
    //Now get rid of the thrown out weights for the particle collection

    // leave these lines in, in case want to move eventually to having no 1-to-1 correspondence between puppi and pf cands
    // if( std::abs(pWeight) < std::numeric_limits<float>::denorm_min() ) continue; // this line seems not to work like it's supposed to...
    // if(std::abs(pWeight) <= 0. ) continue;
  }
  return fWeights;
}
