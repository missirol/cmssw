#ifndef COMMONTOOLS_PUPPI_PUPPICONTAINER_H
#define COMMONTOOLS_PUPPI_PUPPICONTAINER_H

#include <vector>

#include "CommonTools/PileupAlgos/interface/PuppiAlgo.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"

class PuppiContainer {
public:
  PuppiContainer(edm::ParameterSet const& iConfig);
  ~PuppiContainer();
  void initialize(std::vector<PuppiCandidate> const& iRecoObjects);
  void setNPV(int iNPV) { fNPV = iNPV; }

  std::vector<PuppiCandidate> const &pfParticles() const { return fPFParticles; }
  std::vector<PuppiCandidate> const &pvParticles() const { return fChargedPV; }
  std::vector<float> const& puppiWeights();
  const std::vector<float> &puppiRawAlphas() const { return fRawAlphas; }
  const std::vector<float> &puppiAlphas() const { return fVals; }
  const std::vector<float> &puppiAlphasMed() const { return fAlphaMed; }
  const std::vector<float> &puppiAlphasRMS() const { return fAlphaRMS; }

  int puppiNAlgos() const { return fNAlgos; }

protected:
  float goodVar(PuppiCandidate const &iPart, std::vector<PuppiCandidate> const &iParts, int const iOpt, float const iRCone) const;
  void getRMSAvg(int const iOpt,
                 std::vector<PuppiCandidate> const &iConstits,
                 std::vector<PuppiCandidate> const &iParticles,
                 std::vector<PuppiCandidate> const &iChargeParticles);
  void getRawAlphas(int iOpt,
                    std::vector<PuppiCandidate> const &iConstits,
                    std::vector<PuppiCandidate> const &iParticles,
                    std::vector<PuppiCandidate> const &iChargeParticles);
  float getChi2FromdZ(float const iDZ);
  int getPuppiId(float iPt, float iEta);
  float var_within_R(int const iId,
                     std::vector<PuppiCandidate> const& particles,
                     PuppiCandidate const& centre,
                     float const R) const;

  bool fPuppiDiagnostics;
  std::vector<PuppiCandidate> fPFParticles;
  std::vector<PuppiCandidate> fChargedPV;
  std::vector<float> fWeights;
  std::vector<float> fVals;
  std::vector<float> fRawAlphas;
  std::vector<float> fAlphaMed;
  std::vector<float> fAlphaRMS;

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
  int fNPV;
  std::vector<PuppiAlgo> fPuppiAlgo;
};
#endif
