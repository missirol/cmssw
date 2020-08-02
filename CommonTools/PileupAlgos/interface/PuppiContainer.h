#ifndef COMMONTOOLS_PUPPI_PUPPICONTAINER_H
#define COMMONTOOLS_PUPPI_PUPPICONTAINER_H

#include <vector>

#include "CommonTools/PileupAlgos/interface/PuppiAlgo.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"

class PuppiContainer {
public:
  PuppiContainer(edm::ParameterSet const& iConfig);
  ~PuppiContainer();
  void run(std::vector<PuppiCandidate> const& iPuppiCandidates, int const iNPV);

  std::vector<float> const& puppiWeights() const { return fWeights; }
  std::vector<float> const& puppiRawAlphas() const { return fRawAlphas; }
  std::vector<float> const& puppiAlphas() const { return fVals; }
  std::vector<float> const& puppiAlphasMed() const { return fAlphaMed; }
  std::vector<float> const& puppiAlphasRMS() const { return fAlphaRMS; }
  int puppiNAlgos() const { return fNAlgos; }

protected:
  void getRMSAvg(int const iOpt, std::vector<PuppiCandidate> const& iParticles);
  void getRawAlphas(int iOpt, std::vector<PuppiCandidate> const& iParticles);
  float getChi2FromdZ(float const iDZ) const;
  int getPuppiId(float iPt, float iEta);
  float goodVar(PuppiCandidate const &iPart, std::vector<PuppiCandidate> const &iParts, int const iOpt, float const iRCone, bool const iUseCharged) const;
  float var_within_R(int const iId, std::vector<PuppiCandidate> const& particles, PuppiCandidate const& centre, float const R, bool const iUseCharged) const;

  bool fPuppiDiagnostics;
  std::vector<float> fWeights;
  std::vector<float> fVals;
  std::vector<float> fRawAlphas;
  std::vector<float> fAlphaMed;
  std::vector<float> fAlphaRMS;
  std::vector<PuppiAlgo> fPuppiAlgo;

  float fNeutralMinPt;
  float fNeutralSlope;
  float fPuppiWeightCut;
  float fPtMaxPhotons;
  float fEtaMaxPhotons;
  float fPtMaxNeutrals;
  float fPtMaxNeutralsStartSlope;
  int fNAlgos;
  bool fApplyCHS;
  bool fInvert;
  bool fUseExp;
};
#endif
