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

  std::vector<double> const &puppiWeights();
  std::vector<double> const &puppiRawAlphas() const { return fRawAlphas; }
  std::vector<double> const &puppiAlphas() const { return fVals; }
  std::vector<double> const &puppiAlphasMed() const { return fAlphaMed; }
  std::vector<double> const &puppiAlphasRMS() const { return fAlphaRMS; }

  int puppiNAlgos() const { return fNAlgos; }

protected:
  double goodVar(std::vector<PuppiCandidate> const &iCands,
                 uint const iCandIndex0,
                 std::vector<uint> const &iCandIndices,
                 int const iOpt,
                 double const iRCone) const;
  void getRMSAvg(int const iOpt,
                 std::vector<PuppiCandidate> const &iCands,
                 std::vector<uint> const &iCandIndicesSum,
                 std::vector<uint> const &iCandIndicesSumChargedPV);
  void getRawAlphas(int const iOpt,
                    std::vector<PuppiCandidate> const &iCands,
                    std::vector<uint> const &iCandIndicesSum,
                    std::vector<uint> const &iCandIndicesSumChargedPV);
  int getPuppiId(double const iPt, double const iEta);
  double getChi2FromdZ(double const iDZ) const;

  bool fPuppiDiagnostics;
  std::vector<PuppiCandidate> fCands;
  std::vector<uint> fCandIndicesSum;
  std::vector<uint> fCandIndicesSumChargedPV;
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
