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
  int computeWeights(std::vector<RecoObj> const&, int const);

  std::vector<float> const &puppiWeights() const { return fWeights; }
  const std::vector<float> &puppiRawAlphas() const { return fRawAlphas; }
  const std::vector<float> &puppiAlphas() const { return fVals; }
  const std::vector<float> &puppiAlphasMed() const { return fAlphaMed; }
  const std::vector<float> &puppiAlphasRMS() const { return fAlphaRMS; }

  int puppiNAlgos() const { return fNAlgos; }

protected:
  using CandidatesTable = edm::soa::PtEtaRapPhiIdTable;

  float goodVar(CandidatesTable const&, uint const, bool const, int const, float const);
  void getRMSAvg(int iOpt, CandidatesTable const&);
  void getRawAlphas(int iOpt, CandidatesTable const&);

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

#endif
