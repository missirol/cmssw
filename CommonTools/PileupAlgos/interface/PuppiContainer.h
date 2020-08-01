#ifndef COMMONTOOLS_PUPPI_PUPPICONTAINER_H_
#define COMMONTOOLS_PUPPI_PUPPICONTAINER_H_

#include "CommonTools/PileupAlgos/interface/PuppiAlgo.h"
#include "CommonTools/PileupAlgos/interface/RecoObj.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"

class PuppiContainer {
public:
  PuppiContainer(const edm::ParameterSet &iConfig);
  ~PuppiContainer();
  void initialize(const std::vector<RecoObj> &iRecoObjects);
  void setNPV(int iNPV) { fNPV = iNPV; }

  std::vector<PuppiCandidate> const &pfParticles() const { return fPFParticles; }
  std::vector<PuppiCandidate> const &pvParticles() const { return fChargedPV; }
  std::vector<float> const &puppiWeights();
  const std::vector<float> &puppiRawAlphas() { return fRawAlphas; }
  const std::vector<float> &puppiAlphas() { return fVals; }
  // const std::vector<float> puppiAlpha   () {return fAlpha;}
  const std::vector<float> &puppiAlphasMed() { return fAlphaMed; }
  const std::vector<float> &puppiAlphasRMS() { return fAlphaRMS; }

  int puppiNAlgos() { return fNAlgos; }

protected:
  float goodVar(PuppiCandidate const &iPart, std::vector<PuppiCandidate> const &iParts, int iOpt, const float iRCone);
  void getRMSAvg(int iOpt,
                 std::vector<PuppiCandidate> const &iConstits,
                 std::vector<PuppiCandidate> const &iParticles,
                 std::vector<PuppiCandidate> const &iChargeParticles);
  void getRawAlphas(int iOpt,
                    std::vector<PuppiCandidate> const &iConstits,
                    std::vector<PuppiCandidate> const &iParticles,
                    std::vector<PuppiCandidate> const &iChargeParticles);
  float getChi2FromdZ(float const iDZ);
  int getPuppiId(float iPt, float iEta);
  float var_within_R(int iId,
                      const std::vector<PuppiCandidate> &particles,
                      const PuppiCandidate &centre,
                      const float R);

  bool fPuppiDiagnostics;
  const std::vector<RecoObj> *fRecoParticles;
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
  float fPVFrac;
  std::vector<PuppiAlgo> fPuppiAlgo;
};
#endif
