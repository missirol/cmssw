#ifndef CommonTools_PileupAlgos_PuppiContainer_h
#define CommonTools_PileupAlgos_PuppiContainer_h

#include <vector>

#include "CommonTools/PileupAlgos/interface/PuppiAlgo.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"

class PuppiContainer {
public:
  PuppiContainer(edm::ParameterSet const &iConfig);
  ~PuppiContainer();
  void initialize(std::vector<PuppiCandidate> const &iPuppiCandidates);
  void setNPV(int iNPV) { fNPV = iNPV; }

  std::vector<PuppiCandidate> const &pfParticles() const { return fPFParticles; }
  std::vector<PuppiCandidate> const &pvParticles() const { return fChargedPV; }
  std::vector<float> const &puppiWeights();
  std::vector<float> const &puppiRawAlphas() const { return fRawAlphas; }
  std::vector<float> const &puppiAlphas() const { return fVals; }
  std::vector<float> const &puppiAlphasMed() const { return fAlphaMed; }
  std::vector<float> const &puppiAlphasRMS() const { return fAlphaRMS; }

  int puppiNAlgos() const { return fNAlgos; }

protected:
  float goodVar(PuppiCandidate const &iPart,
                 std::vector<PuppiCandidate> const &iParticles,
                 int const iOpt,
                 float const iRCone) const;
  void getRMSAvg(int const iOpt,
                 std::vector<PuppiCandidate> const &iParticles,
                 std::vector<PuppiCandidate> const &iParticlesXXX,
                 std::vector<PuppiCandidate> const &iChargeParticles);
  void getRawAlphas(int const iOpt,
                    std::vector<PuppiCandidate> const &iParticles,
                    std::vector<PuppiCandidate> const &iParticlesXXX,
                    std::vector<PuppiCandidate> const &iChargeParticles);
  int getPuppiId(float const iPt, float const iEta);
  float getChi2FromdZ(float const iDZ) const;

  bool fPuppiDiagnostics;
  std::vector<PuppiCandidate> fPFParticles;
  std::vector<PuppiCandidate> fPFParticlesXXX;
  std::vector<PuppiCandidate> fChargedPV;
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
  int fNPV;
  bool fApplyCHS;
  bool fInvert;
  bool fUseExp;
};
#endif
