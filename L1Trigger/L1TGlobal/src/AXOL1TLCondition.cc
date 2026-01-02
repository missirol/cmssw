/**
 * \class AXOL1TLCondition
 *
 *
 * Description: evaluation of a condition for axol1tl anomaly detection algorithm
 *
 * Author: Melissa Quinnan
 *
 **/
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include <utility>

#include "ap_fixed.h"

#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "L1Trigger/L1TGlobal/interface/AXOL1TLCondition.h"
#include "L1Trigger/L1TGlobal/interface/AXOL1TLTemplate.h"
#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"
#include "L1Trigger/L1TGlobal/interface/GlobalBoard.h"

namespace {
  // template function for reading results
  template <typename InputType, typename ResultType, typename LossType>
  LossType readResult(hls4mlEmulator::ModelWrapper const& modelWrapper, InputType inputs[]) {
    // model outputs a pair of the (result vector, loss)
    std::pair<ResultType, LossType> ADModelResult;
    modelWrapper.run_inference(inputs, &ADModelResult);
    return ADModelResult.second;
  }
}  // namespace

l1t::AXOL1TLCondition::AXOL1TLCondition()
    : ConditionEvaluation(), m_gtAXOL1TLTemplate{nullptr}, m_gtGTB{nullptr}, m_model_wrapper{}, m_saved_score{0} {}

l1t::AXOL1TLCondition::AXOL1TLCondition(const GlobalCondition* axol1tlTemplate, const GlobalBoard* ptrGTB)
    : ConditionEvaluation(),
      m_gtAXOL1TLTemplate(static_cast<const AXOL1TLTemplate*>(axol1tlTemplate)),
      m_gtGTB(ptrGTB),
      m_model_wrapper{kModelNamePrefix + m_gtAXOL1TLTemplate->modelVersion()},
      m_saved_score{0} {}

// copy constructor
void l1t::AXOL1TLCondition::copy(const l1t::AXOL1TLCondition& cp) {
  m_condMaxNumberObjects = cp.condMaxNumberObjects();
  m_condLastResult = cp.condLastResult();
  m_combinationsInCond = cp.getCombinationsInCond();
  m_verbosity = cp.m_verbosity;

  m_gtAXOL1TLTemplate = cp.gtAXOL1TLTemplate();
  m_gtGTB = cp.gtGTB();
  m_model_wrapper.reset(cp.model_name());
  m_saved_score = cp.getScore();
}

l1t::AXOL1TLCondition::AXOL1TLCondition(const l1t::AXOL1TLCondition& cp) : ConditionEvaluation() { copy(cp); }

l1t::AXOL1TLCondition& l1t::AXOL1TLCondition::operator=(const l1t::AXOL1TLCondition& cp) {
  copy(cp);
  return *this;
}

const bool l1t::AXOL1TLCondition::evaluateCondition(const int bxEval) const {
  int const useBx = bxEval + m_gtAXOL1TLTemplate->condRelativeBx();

  // pointers to objects
  const BXVector<const l1t::EtSum*>* candEtSumVec = m_gtGTB->getCandL1EtSum();
  const BXVector<const l1t::L1Candidate*>* candEGVec = m_gtGTB->getCandL1EG();
  const BXVector<const l1t::Muon*>* candMuVec = m_gtGTB->getCandL1Mu();
  const BXVector<const l1t::L1Candidate*>* candJetVec = m_gtGTB->getCandL1Jet();

  int const NEtSums = 1;
  int const NEgammas = 4;
  int const NMuons = 4;
  int const NJets = 10;

  // number of input features: #objects * 3 (for et, eta, phi)
  // total: (1 + 4 + 4 + 10) * 3 = 57
  int const EtSumVecSize = 3 * NEtSums;
  int const EGVecSize = 3 * NEgammas;
  int const MuVecSize = 3 * NMuons;
  int const JVecSize = 3 * NJets;

  int const NInputs = EtSumVecSize + EGVecSize + MuVecSize + JVecSize;

  // types of inputs and outputs
  typedef ap_fixed<18, 13> inputtype;
  typedef ap_ufixed<18, 14> losstype;

  // arrays of input features per object
  inputtype EtSumInput[EtSumVecSize];
  inputtype EgammaInput[EGVecSize];
  inputtype MuInput[MuVecSize];
  inputtype JetInput[JVecSize];
  inputtype ADModelInput[NInputs] = {};

  // output object
  losstype loss;

  // check number of input objects we actually have (muons, jets etc)
  int const NCandEtSum = candEtSumVec->size(useBx);
  int const NCandEG = candEGVec->size(useBx);
  int const NCandMu = candMuVec->size(useBx);
  int const NCandJet = candJetVec->size(useBx);

  // initialize arrays to zero (std::fill(first, last, value);)
  inputtype const fillzero = 0.0;
  std::fill(EtSumInput, EtSumInput + EtSumVecSize, fillzero);
  std::fill(EgammaInput, EgammaInput + EGVecSize, fillzero);
  std::fill(MuInput, MuInput + MuVecSize, fillzero);
  std::fill(JetInput, JetInput + JVecSize, fillzero);
  std::fill(ADModelInput, ADModelInput + NInputs, fillzero);

  // then fill the object arrays
  // NOTE assume candidates are already sorted by pt

  // loop over EtSums first
  for (int iEtSum = 0; iEtSum < NCandEtSum; iEtSum++) {
    if (iEtSum < NEtSums and candEtSumVec->at(useBx, iEtSum)->getType() == l1t::EtSum::EtSumType::kMissingEt) {
      // have to do hwPt/2 in order to match original et inputs
      EtSumInput[0 + (3 * iEtSum)] = candEtSumVec->at(useBx, iEtSum)->hwPt() / 2;
      // leave EtSumInput[1 + (3 * iEtSum)] (eta) at zero
      EtSumInput[2 + (3 * iEtSum)] = candEtSumVec->at(useBx, iEtSum)->hwPhi();
    }
  }

  // next egammas
  for (int iEG = 0; iEG < NCandEG; iEG++) {
    if (iEG < NEgammas) {
      // have to do hwPt/2 in order to match original et inputs
      EgammaInput[0 + (3 * iEG)] = candEGVec->at(useBx, iEG)->hwPt() / 2;
      EgammaInput[1 + (3 * iEG)] = candEGVec->at(useBx, iEG)->hwEta();
      EgammaInput[2 + (3 * iEG)] = candEGVec->at(useBx, iEG)->hwPhi();
    }
  }

  // next muons
  for (int iMu = 0; iMu < NCandMu; iMu++) {
    if (iMu < NMuons) {
      // have to do hwPt/2 in order to match original et inputs
      MuInput[0 + (3 * iMu)] = candMuVec->at(useBx, iMu)->hwPt() / 2;
      MuInput[1 + (3 * iMu)] = candMuVec->at(useBx, iMu)->hwEtaAtVtx();
      MuInput[2 + (3 * iMu)] = candMuVec->at(useBx, iMu)->hwPhiAtVtx();
    }
  }

  // next jets
  for (int iJet = 0; iJet < NCandJet; iJet++) {
    if (iJet < NJets) {
      // have to do hwPt/2 in order to match original et inputs
      JetInput[0 + (3 * iJet)] = candJetVec->at(useBx, iJet)->hwPt() / 2;
      JetInput[1 + (3 * iJet)] = candJetVec->at(useBx, iJet)->hwEta();
      JetInput[2 + (3 * iJet)] = candJetVec->at(useBx, iJet)->hwPhi();
    }
  }

  // now put it all together-> EtSum+EGamma+Muon+Jet into ADModelInput
  int index = 0;
  for (int idET = 0; idET < EtSumVecSize; idET++) {
    ADModelInput[index++] = EtSumInput[idET];
  }
  for (int idEG = 0; idEG < EGVecSize; idEG++) {
    ADModelInput[index++] = EgammaInput[idEG];
  }
  for (int idMu = 0; idMu < MuVecSize; idMu++) {
    ADModelInput[index++] = MuInput[idMu];
  }
  for (int idJ = 0; idJ < JVecSize; idJ++) {
    ADModelInput[index++] = JetInput[idJ];
  }

  // now run the inference
  if (m_model_wrapper.model_name() == "GTADModel_v3" or m_model_wrapper.model_name() == "GTADModel_v4") {
    using resulttype = std::array<ap_fixed<10, 7, AP_RND_CONV, AP_SAT>, 8>;
    loss = readResult<inputtype, resulttype, losstype>(m_model_wrapper, ADModelInput);
  } else {
    using resulttype = ap_fixed<18, 14, AP_RND_CONV, AP_SAT>;
    loss = readResult<inputtype, resulttype, losstype>(m_model_wrapper, ADModelInput);
  }

  // scaling to match threshold
  float const score = loss.to_float() * 16;

  // save score to class variable in case score saving needed
  setScore(score);

  // number of objects/thresholds to check
  int const nObjInCond = m_gtAXOL1TLTemplate->nrObjects();

  // number of conditions: there is only one
  int const iCondition = 0;

  if (iCondition >= nObjInCond || iCondition < 0) {
    return false;
  }

  AXOL1TLTemplate::ObjectParameter const objPar = (*(m_gtAXOL1TLTemplate->objectParameter()))[iCondition];

  // condGEqVal indicates the operator used for the condition (>=, =): true for >=
  bool const condGEqVal = m_gtAXOL1TLTemplate->condGEq();

  bool const condResult = checkCut(objPar.minAXOL1TLThreshold, score, condGEqVal);

  // return result
  return condResult;
}

void l1t::AXOL1TLCondition::print(std::ostream& myCout) const {
  myCout << "Dummy Print for AXOL1TLCondition" << std::endl;
  m_gtAXOL1TLTemplate->print(myCout);

  ConditionEvaluation::print(myCout);
}
