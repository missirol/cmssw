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

  std::vector<PuppiCandidate> const &pfParticles() const { return fPFParticles; }
  std::vector<PuppiCandidate> const &pvParticles() const { return fChargedPV; }
  std::vector<double> const &puppiWeights(const int);
  const std::vector<double> &puppiRawAlphas() { return fRawAlphas; }
  const std::vector<double> &puppiAlphas() { return fVals; }
  // const std::vector<double> puppiAlpha   () {return fAlpha;}
  const std::vector<double> &puppiAlphasMed() { return fAlphaMed; }
  const std::vector<double> &puppiAlphasRMS() { return fAlphaRMS; }

  int puppiNAlgos() { return fNAlgos; }

protected:
  double cached_dr2(uint, uint) const;
  double goodVar(const uint, std::vector<PuppiCandidate> const &iParts, int iOpt, const double iRCone);
  void getRMSAvg(int iOpt,
                 std::vector<PuppiCandidate> const &iConstits,
                 std::vector<PuppiCandidate> const &iParticles,
                 std::vector<PuppiCandidate> const &iChargeParticles);
  void getRawAlphas(int iOpt,
                    std::vector<PuppiCandidate> const &iConstits,
                    std::vector<PuppiCandidate> const &iParticles,
                    std::vector<PuppiCandidate> const &iChargeParticles);
  double getChi2FromdZ(double iDZ);
  int getPuppiId(float iPt, float iEta);
  double var_within_R(int iId, const std::vector<PuppiCandidate> &particles, const uint, const double R);

  bool fPuppiDiagnostics;
  const std::vector<RecoObj> *fRecoParticles;
  std::vector<PuppiCandidate> fPFParticles;
  std::vector<PuppiCandidate> fChargedPV;
  std::vector<double> fWeights;
  std::vector<double> fVals;
  std::vector<double> fRawAlphas;
  std::vector<double> fAlphaMed;
  std::vector<double> fAlphaRMS;
  std::vector<double> fDeltaR2;

  bool fApplyCHS;
  bool fInvert;
  bool fUseExp;
  double fNeutralMinPt;
  double fNeutralSlope;
  double fPuppiWeightCut;
  double fPtMaxPhotons;
  double fEtaMaxPhotons;
  double fPtMaxNeutrals;
  double fPtMaxNeutralsStartSlope;
  int fNAlgos;
  std::vector<PuppiAlgo> fPuppiAlgo;
};
#endif
