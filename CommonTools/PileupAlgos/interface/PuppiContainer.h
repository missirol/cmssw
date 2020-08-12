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
  void setNPV(int const iNPV) { fNPV = iNPV; }

  std::vector<PuppiCandidate> const &puppiCandidates() const { return fCands; }
  std::vector<float> const &puppiWeights();
  std::vector<float> const &puppiRawAlphas() const { return fRawAlphas; }
  std::vector<float> const &puppiAlphas() const { return fVals; }
  std::vector<float> const &puppiAlphasMed() const { return fAlphaMed; }
  std::vector<float> const &puppiAlphasRMS() const { return fAlphaRMS; }

  int puppiNAlgos() const { return fNAlgos; }

protected:
  float goodVar(PuppiCandidate const &iPuppiCand_0,
                std::vector<PuppiCandidate> const &iPuppiCands,
                int const iOpt,
                float const iRCone) const;
  void getRMSAvg(int const iOpt,
                 std::vector<PuppiCandidate> const &iPuppiCands,
                 std::vector<PuppiCandidate> const &iPuppiCandsForVar,
                 std::vector<PuppiCandidate> const &iPuppiCandsForVarChargedPV);
  void getRawAlphas(int const iOpt,
                    std::vector<PuppiCandidate> const &iPuppiCands,
                    std::vector<PuppiCandidate> const &iPuppiCandsForVar,
                    std::vector<PuppiCandidate> const &iPuppiCandsForVarChargedPV);
  int getPuppiId(float const iPt, float const iEta);
  float getChi2FromdZ(float const iDZ) const;

  bool fPuppiDiagnostics;
  std::vector<PuppiCandidate> fCands;
  std::vector<PuppiCandidate> fCandsForVar;
  std::vector<PuppiCandidate> fCandsForVarChargedPV;
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
