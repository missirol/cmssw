#ifndef CommonTools_PileupAlgos_PuppiAlgo_h
#define CommonTools_PileupAlgos_PuppiAlgo_h

#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"

class PuppiAlgo {
public:
  PuppiAlgo(edm::ParameterSet &iConfig);
  ~PuppiAlgo();
  static void fillDescriptionsPuppiAlgo(edm::ParameterSetDescription &desc);
  //Computing Mean and RMS
  void reset();
  void fixAlgoEtaBin(int i_eta);
  void add(PuppiCandidate const &iParticle, double const iVal, uint const iAlgo);
  void computeMedRMS(uint const iAlgo);
  //Get the Weight
  double compute(std::vector<double> const &iVals, double const iChi2) const;
  std::vector<double> const &alphas() const { return fPups; }
  //Helpers
  inline int etaBins() const { return fEtaMin.size(); }
  inline double etaMin(int i) const { return fEtaMin[i]; }
  inline double etaMax(int i) const { return fEtaMax[i]; }
  inline double ptMin() const { return cur_PtMin; }

  inline int numAlgos() const { return fNAlgos; }
  inline int algoId(uint iAlgo) const { return fAlgoId.at(iAlgo); }
  inline bool isCharged(uint iAlgo) const { return fCharged.at(iAlgo); }
  inline double coneSize(uint iAlgo) const { return fConeSize.at(iAlgo); }
  inline double neutralPt(int iNPV) const { return cur_NeutralPtMin + iNPV * cur_NeutralPtSlope; }

  inline double rms() const { return cur_RMS; }
  inline double median() const { return cur_Med; }

  inline double etaMaxExtrap() const { return fEtaMaxExtrap; }

private:
  uint fNAlgos;
  std::vector<double> fEtaMax;
  std::vector<double> fEtaMin;
  std::vector<double> fPtMin;
  std::vector<double> fNeutralPtMin;
  std::vector<double> fNeutralPtSlope;

  std::vector<double> fRMSEtaSF;
  std::vector<double> fMedEtaSF;
  double fEtaMaxExtrap;

  double cur_PtMin;
  double cur_NeutralPtMin;
  double cur_NeutralPtSlope;
  double cur_RMS;
  double cur_Med;

  std::vector<double> fRMS;                          // raw RMS per algo
  std::vector<double> fMedian;                       // raw Median per algo
  std::vector<std::vector<double> > fRMS_perEta;     // final RMS used after eta corrections
  std::vector<std::vector<double> > fMedian_perEta;  // final Med used after eta corrections

  std::vector<double> fPups;
  std::vector<double> fPupsPV;
  std::vector<int> fAlgoId;
  std::vector<bool> fCharged;
  std::vector<bool> fAdjust;
  std::vector<int> fCombId;
  std::vector<double> fConeSize;
  std::vector<double> fRMSPtMin;
  std::vector<double> fRMSScaleFactor;
  std::vector<double> fMean;
  std::vector<int> fNCount;
};

#endif
