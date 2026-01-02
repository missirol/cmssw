/**
 * \class TOPOCondition
 *
 *
 * Description: evaluation of a condition for TOPO anomaly detection algorithm
 *
 * Author: Melissa Quinnan, Lukas Ebeling, Artur Lobanov
 *
 **/
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>

#include "ap_fixed.h"

#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"
#include "L1Trigger/L1TGlobal/interface/GlobalBoard.h"
#include "L1Trigger/L1TGlobal/interface/TOPOCondition.h"
#include "L1Trigger/L1TGlobal/interface/TOPOTemplate.h"

l1t::TOPOCondition::TOPOCondition()
    : ConditionEvaluation(), m_gtTOPOTemplate{nullptr}, m_gtGTB{nullptr}, m_model_wrapper{} {}

l1t::TOPOCondition::TOPOCondition(const GlobalCondition* topoTemplate, const GlobalBoard* ptrGTB)
    : ConditionEvaluation(),
      m_gtTOPOTemplate(static_cast<const TOPOTemplate*>(topoTemplate)),
      m_gtGTB(ptrGTB),
      m_model_wrapper{kModelNamePrefix + m_gtTOPOTemplate->modelVersion()} {}

void l1t::TOPOCondition::copy(const l1t::TOPOCondition& cp) {
  m_condMaxNumberObjects = cp.condMaxNumberObjects();
  m_condLastResult = cp.condLastResult();
  m_combinationsInCond = cp.getCombinationsInCond();
  m_verbosity = cp.m_verbosity;

  m_gtTOPOTemplate = cp.gtTOPOTemplate();
  m_gtGTB = cp.gtGTB();
  m_model_wrapper.reset(cp.model_name());
}

l1t::TOPOCondition::TOPOCondition(const l1t::TOPOCondition& cp) : ConditionEvaluation() { copy(cp); }

l1t::TOPOCondition& l1t::TOPOCondition::operator=(const l1t::TOPOCondition& cp) {
  copy(cp);
  return *this;
}

const bool l1t::TOPOCondition::evaluateCondition(const int bxEval) const {
  int const useBx = bxEval + m_gtTOPOTemplate->condRelativeBx();

  // pointers to objects
  const BXVector<const l1t::EtSum*>* candEtSumVec = m_gtGTB->getCandL1EtSum();
  const BXVector<const l1t::L1Candidate*>* candEGVec = m_gtGTB->getCandL1EG();
  const BXVector<const l1t::Muon*>* candMuVec = m_gtGTB->getCandL1Mu();
  const BXVector<const l1t::L1Candidate*>* candJetVec = m_gtGTB->getCandL1Jet();

  // number of objects and input features
  // total: (1*1 + 0*3 + 2*4 + 4*3) = 21
  int const NEtSums = 1;
  int const NEgammas = 0;
  int const NMuons = 2;
  int const NJets = 4;

  int const EtSumVecSize = NEtSums * 1;
  int const EGVecSize = NEgammas * 3;
  int const MuVecSize = NMuons * 4;
  int const JVecSize = NJets * 3;

  int const NInputs = EtSumVecSize + EGVecSize + MuVecSize + JVecSize;

  // types of inputs and outputs
  typedef ap_fixed<23, 23> inputtype;
  typedef ap_fixed<16, 6> losstype;

  // arrays of input features per object
  double EtSumInput[EtSumVecSize];
  double EgammaInput[EGVecSize];
  double MuInput[MuVecSize];
  double JetInput[JVecSize];
  inputtype ModelInput[NInputs];

  // output object
  losstype loss;

  // check number of input objects we actually have (muons, jets etc)
  int const NCandEtSum = candEtSumVec->size(useBx);
  int const NCandEG = candEGVec->size(useBx);
  int const NCandMu = candMuVec->size(useBx);
  int const NCandJet = candJetVec->size(useBx);

  // initialize arrays to zero (std::fill(first, last, value))
  inputtype const fillzero = 0.0;
  std::fill(EtSumInput, EtSumInput + EtSumVecSize, fillzero);
  std::fill(EgammaInput, EgammaInput + EGVecSize, fillzero);
  std::fill(MuInput, MuInput + MuVecSize, fillzero);
  std::fill(JetInput, JetInput + JVecSize, fillzero);
  std::fill(ModelInput, ModelInput + NInputs, fillzero);

  // fill EtSum array
  for (int iEtSum = 0; iEtSum < NCandEtSum; iEtSum++) {
    if (iEtSum < NEtSums and candEtSumVec->at(useBx, iEtSum)->getType() == l1t::EtSum::EtSumType::kTotalHt) {
      EtSumInput[0 + (1 * iEtSum)] = candEtSumVec->at(useBx, iEtSum)->hwPt();
    }
  }

  // next egammas
  for (int iEG = 0; iEG < NCandEG; iEG++) {
    if (iEG < NEgammas) {
      EgammaInput[0 + (3 * iEG)] = candEGVec->at(useBx, iEG)->hwPt();
      EgammaInput[1 + (3 * iEG)] = candEGVec->at(useBx, iEG)->hwEta();
      EgammaInput[2 + (3 * iEG)] = candEGVec->at(useBx, iEG)->hwPhi();
    }
  }

  // next muons
  for (int iMu = 0; iMu < NCandMu; iMu++) {
    if (iMu < NMuons) {
      MuInput[0 + (4 * iMu)] = candMuVec->at(useBx, iMu)->hwPt();
      MuInput[1 + (4 * iMu)] = candMuVec->at(useBx, iMu)->hwEtaAtVtx();
      MuInput[2 + (4 * iMu)] = candMuVec->at(useBx, iMu)->hwPhiAtVtx();
      MuInput[3 + (4 * iMu)] = candMuVec->at(useBx, iMu)->hwQual();
    }
  }

  // next jets
  for (int iJet = 0; iJet < NCandJet; iJet++) {
    if (iJet < NJets) {
      JetInput[0 + (3 * iJet)] = candJetVec->at(useBx, iJet)->hwPt();
      JetInput[1 + (3 * iJet)] = candJetVec->at(useBx, iJet)->hwEta();
      JetInput[2 + (3 * iJet)] = candJetVec->at(useBx, iJet)->hwPhi();
    }
  }

  // now put it all together -> EtSum+EGamma+Muon+Jet into ModelInput
  int index = 0;
  for (int idET = 0; idET < EtSumVecSize; idET++) {
    ModelInput[index++] = EtSumInput[idET];
  }
  for (int idEG = 0; idEG < EGVecSize; idEG++) {
    ModelInput[index++] = EgammaInput[idEG];
  }
  for (int idMu = 0; idMu < MuVecSize; idMu++) {
    ModelInput[index++] = MuInput[idMu];
  }
  for (int idJ = 0; idJ < JVecSize; idJ++) {
    ModelInput[index++] = JetInput[idJ];
  }

  // now run the inference
  m_model_wrapper.run_inference(ModelInput, &loss);

  float const score = loss.to_float() * 1023;

  // number of objects/thresholds to check
  int const nObjInCond = m_gtTOPOTemplate->nrObjects();

  // number of conditions: there is only one
  int const iCondition = 0;

  if (iCondition >= nObjInCond || iCondition < 0) {
    return false;
  }

  TOPOTemplate::ObjectParameter const objPar = (*(m_gtTOPOTemplate->objectParameter()))[iCondition];

  // condGEqVal indicates the operator used for the condition (>=, =): true for >=
  bool const condGEqVal = m_gtTOPOTemplate->condGEq();

  bool const condResult = checkCut(objPar.minTOPOThreshold, score, condGEqVal);

  // return result
  return condResult;
}

void l1t::TOPOCondition::print(std::ostream& myCout) const {
  myCout << "Dummy Print for TOPOCondition" << std::endl;
  m_gtTOPOTemplate->print(myCout);

  ConditionEvaluation::print(myCout);
}
