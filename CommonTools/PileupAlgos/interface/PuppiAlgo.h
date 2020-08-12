#ifndef CommonTools_PileupAlgos_PuppiAlgo_h
#define CommonTools_PileupAlgos_PuppiAlgo_h

#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"

class PuppiAlgo {
public:
  PuppiAlgo(edm::ParameterSet const &iConfig);
  ~PuppiAlgo();
  static void fillDescriptionsPuppiAlgo(edm::ParameterSetDescription &desc);
  //Computing Mean and RMS
  void reset();
  void fixAlgoEtaBin(int i_eta);
  void add(PuppiCandidate const &iParticle, float const iVal, uint const iAlgo);
  void computeMedRMS(uint const iAlgo);
  //Get the Weight
  float compute(std::vector<float> const &iVals, float const iChi2) const;
  std::vector<float> const &alphas() const { return fPups; }
  //Helpers
  inline int etaBins() const { return fEtaMin.size(); }
  inline float etaMin(int i) const { return fEtaMin[i]; }
  inline float etaMax(int i) const { return fEtaMax[i]; }
  inline float ptMin() const { return cur_PtMin; }

  inline int numAlgos() const { return fNAlgos; }
  inline int algoId(uint iAlgo) const { return fAlgoId.at(iAlgo); }
  inline bool isCharged(uint iAlgo) const { return fCharged.at(iAlgo); }
  inline float coneSize(uint iAlgo) const { return fConeSize.at(iAlgo); }
  inline float neutralPt(int iNPV) const { return cur_NeutralPtMin + iNPV * cur_NeutralPtSlope; }

  inline float rms() const { return cur_RMS; }
  inline float median() const { return cur_Med; }

  inline float etaMaxExtrap() const { return fEtaMaxExtrap; }

private:
  uint fNAlgos;
  std::vector<double> fEtaMax;
  std::vector<double> fEtaMin;
  std::vector<double> fPtMin;
  std::vector<double> fNeutralPtMin;
  std::vector<double> fNeutralPtSlope;
  std::vector<double> fRMSEtaSF;
  std::vector<double> fMedEtaSF;
  float fEtaMaxExtrap;

  float cur_PtMin;
  float cur_NeutralPtMin;
  float cur_NeutralPtSlope;
  float cur_RMS;
  float cur_Med;

  std::vector<float> fRMS;                          // raw RMS per algo
  std::vector<float> fMedian;                       // raw Median per algo
  std::vector<std::vector<float> > fRMS_perEta;     // final RMS used after eta corrections
  std::vector<std::vector<float> > fMedian_perEta;  // final Med used after eta corrections

  std::vector<float> fPups;
  std::vector<float> fPupsPV;
  std::vector<int> fAlgoId;
  std::vector<bool> fCharged;
  std::vector<bool> fAdjust;
  std::vector<int> fCombId;
  std::vector<float> fConeSize;
  std::vector<float> fRMSPtMin;
  std::vector<float> fRMSScaleFactor;
  std::vector<float> fMean;
  std::vector<int> fNCount;
};

#endif
