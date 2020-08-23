#ifndef COMMONTOOLS_PUPPI_PUPPICONTAINER_H_
#define COMMONTOOLS_PUPPI_PUPPICONTAINER_H_

#include "CommonTools/PileupAlgos/interface/PuppiAlgo.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidates.h"
#include "CommonTools/PileupAlgos/interface/RecoObj.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/ProbFunc.h"
#include "TMath.h"
#include <iostream>
#include <cmath>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"


class PuppiContainer {
public:
  PuppiContainer(const edm::ParameterSet &iConfig);
  ~PuppiContainer();
  int fillValues(std::vector<RecoObj> const&, int const);

  std::vector<float> const &puppiWeights() const { return fWeights; }
  const std::vector<float> &puppiRawAlphas() { return fRawAlphas; }
  const std::vector<float> &puppiAlphas() { return fVals; }
  const std::vector<float> &puppiAlphasMed() { return fAlphaMed; }
  const std::vector<float> &puppiAlphasRMS() { return fAlphaRMS; }

  int puppiNAlgos() { return fNAlgos; }

protected:
  template<unsigned int N>
  float goodVar(PuppiCandidates<N> const&, uint const, bool const, int const, float const);

  template<unsigned int N>
  void getRMSAvg(int iOpt, PuppiCandidates<N> const&);

  template<unsigned int N>
  void getRawAlphas(int iOpt, PuppiCandidates<N> const&);

  float getChi2FromdZ(float iDZ);
  int getPuppiId(float iPt, float iEta);

  std::vector<float> fWeights;
  std::vector<float> fVals;
  std::vector<float> fRawAlphas;
  std::vector<float> fAlphaMed;
  std::vector<float> fAlphaRMS;

  bool fPuppiDiagnostics;
  bool fApplyCHS;
  bool fInvert;
  bool fUseExp;
  float fNeutralMinPt;
  float fNeutralSlope;
  float fPuppiWeightCut;
  float fPtMaxPhotons;
  float fEtaMaxPhotons;
  float fPtMaxNeutrals;
  float fPtMaxNeutralsStartSlope;
  int fNAlgos;
  std::vector<PuppiAlgo> fPuppiAlgo;
};

template<unsigned int N>
float PuppiContainer::goodVar(PuppiCandidates<N> const &candidates,
                              uint const candIdx,
                              bool const useLVCands,
                              int const iId,
                              float const R) {
  if (iId == -1)
    return 1.f;

  float const r2 = R * R;
  float var = 0.f;

  for (uint idx=0; idx<candidates.size(); ++idx) {
    if(idx == candIdx or candidates.id[candIdx] == 3 or (useLVCands and candidates.id[candIdx] == 1))
      continue;

    if (std::abs(candidates.rapidity[idx] - candidates.rapidity[candIdx]) < R) {
      auto const dr2y = reco::deltaR2(candidates.rapidity[idx], candidates.phi[idx], candidates.rapidity[candIdx], candidates.phi[candIdx]);
      if (dr2y < r2) {
        auto const dr2 = reco::deltaR2(candidates.eta[idx], candidates.phi[idx], candidates.eta[candIdx], candidates.phi[candIdx]);
        auto const pt = candidates.pt[idx];
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

  if ((var != 0.) and ((iId == 0) or (iId == 3) or (iId == 5)))
    var = log(var);
  else if (iId == 1)
    var += candidates.pt[candIdx];

  return var;
}

template<unsigned int N>
void PuppiContainer::getRMSAvg(int iOpt, PuppiCandidates<N> const &iConstits) {
  for (unsigned int i0 = 0; i0 < iConstits.size(); i0++) {
    //Calculate the Puppi Algo to use                                                                                                                                                                      
    int pPupId = getPuppiId(iConstits[i0].pt, iConstits[i0].eta);
    if (pPupId == -1 || fPuppiAlgo[pPupId].numAlgos() <= iOpt) {
      fVals.push_back(-1);
      continue;
    }
    //Get the Puppi Sub Algo (given iteration)                                                                                                                                                             
    int pAlgo = fPuppiAlgo[pPupId].algoId(iOpt);
    bool pCharged = fPuppiAlgo[pPupId].isCharged(iOpt);
    float pCone = fPuppiAlgo[pPupId].coneSize(iOpt);
    // compute the Puppi metric:                                                                                                                                                                           
    //  - calculate goodVar only for candidates that (1) will not be assigned a predefined weight (e.g 0, 1),                                                                                              
    //    or (2) are required for computations inside puppi-algos (see call to PuppiAlgo::add below)                                                                                                       
    float pVal = -1;
    bool const getsDefaultWgtIfApplyCHS = iConstits[i0].id == 1 or iConstits[i0].id == 2;
    if (not((fApplyCHS and getsDefaultWgtIfApplyCHS) or iConstits[i0].id == 3) or
        (std::abs(iConstits[i0].eta) < fPuppiAlgo[pPupId].etaMaxExtrap() and getsDefaultWgtIfApplyCHS)) {
      pVal = goodVar(iConstits, i0, pCharged, pAlgo, pCone);
    }
    fVals.push_back(pVal);

    if (!edm::isFinite(pVal)) {
      LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iConstits[i0].pt << " -- " << iConstits[i0].eta;
      continue;
    }

    // code added by Nhan: now instead for every algorithm give it all the particles                                                                                                                       
    for (int i1 = 0; i1 < fNAlgos; i1++) {
      // skip cands outside of algo's etaMaxExtrap, as they would anyway be ignored inside PuppiAlgo::add (see end of the block)                                                                           
      if (not(std::abs(iConstits[i0].eta) < fPuppiAlgo[i1].etaMaxExtrap() and getsDefaultWgtIfApplyCHS))
        continue;

      auto curVal = pVal;
      // recompute goodVar if algo has changed                                                                                                                                                             
      if (i1 != pPupId) {
        pAlgo = fPuppiAlgo[i1].algoId(iOpt);
        pCharged = fPuppiAlgo[i1].isCharged(iOpt);
        pCone = fPuppiAlgo[i1].coneSize(iOpt);
        curVal = goodVar(iConstits, i0, pCharged, pAlgo, pCone);
      }

      fPuppiAlgo[i1].add(curVal, iOpt, iConstits.pt[i0], iConstits.eta[i0], iConstits.id[i0]);
    }
  }

  for (int i0 = 0; i0 < fNAlgos; i0++)
    fPuppiAlgo[i0].computeMedRMS(iOpt);
}

template<unsigned int N>
void PuppiContainer::getRawAlphas(int iOpt, PuppiCandidates<N> const &iConstits) {
  for (int j0 = 0; j0 < fNAlgos; j0++) {
    for (unsigned int i0 = 0; i0 < iConstits.size(); i0++) {
      //Get the Puppi Sub Algo (given iteration)                                                                                                                                                           
      int pAlgo = fPuppiAlgo[j0].algoId(iOpt);
      bool pCharged = fPuppiAlgo[j0].isCharged(iOpt);
      float pCone = fPuppiAlgo[j0].coneSize(iOpt);
      //Compute the Puppi Metric                                                                                                                                                                           
      float const pVal = goodVar(iConstits, i0, pCharged, pAlgo, pCone);
      fRawAlphas.push_back(pVal);
      if (!edm::isFinite(pVal)) {
        LogDebug("NotFound") << "====> Value is Nan " << pVal << " == " << iConstits[i0].pt << " -- "
                             << iConstits[i0].eta;
        continue;
      }
    }
  }
}

#endif
