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
  std::vector<double> const &puppiWeights();
  std::vector<double> const &puppiRawAlphas() const { return fRawAlphas; }
  std::vector<double> const &puppiAlphas() const { return fVals; }
  std::vector<double> const &puppiAlphasMed() const { return fAlphaMed; }
  std::vector<double> const &puppiAlphasRMS() const { return fAlphaRMS; }

  int puppiNAlgos() const { return fNAlgos; }

protected:
  double goodVar(PuppiCandidate const &iPart,
                 std::vector<PuppiCandidate> const &iParticles,
                 int const iOpt,
                 double const iRCone) const;
  void getRMSAvg(int const iOpt,
                 std::vector<PuppiCandidate> const &iParticles,
                 std::vector<PuppiCandidate> const &iChargeParticles);
  void getRawAlphas(int const iOpt,
                    std::vector<PuppiCandidate> const &iParticles,
                    std::vector<PuppiCandidate> const &iChargeParticles);
  int getPuppiId(double const iPt, double const iEta);
  double getChi2FromdZ(double const iDZ) const;

  bool fPuppiDiagnostics;
  std::vector<PuppiCandidate> fPFParticles;
  std::vector<PuppiCandidate> fChargedPV;
  std::vector<double> fWeights;
  std::vector<double> fVals;
  std::vector<double> fRawAlphas;
  std::vector<double> fAlphaMed;
  std::vector<double> fAlphaRMS;
  std::vector<PuppiAlgo> fPuppiAlgo;

  double fNeutralMinPt;
  double fNeutralSlope;
  double fPuppiWeightCut;
  double fPtMaxPhotons;
  double fEtaMaxPhotons;
  double fPtMaxNeutrals;
  double fPtMaxNeutralsStartSlope;
  int fNAlgos;
  int fNPV;
  bool fApplyCHS;
  bool fInvert;
  bool fUseExp;
};
#endif
