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
  //template function for reading results
  template <typename InputType, typename ResultType, typename LossType>
  LossType readResult(hls4mlEmulator::ModelWrapper const& modelWrapper, InputType const inputs[]) {
    // model outputs a pair of the (result vector, loss)
    std::pair<ResultType, LossType> ADModelResult;
    modelWrapper.run_inference(inputs, &ADModelResult);
    return ADModelResult.second;
  }
}  // namespace

l1t::AXOL1TLCondition::AXOL1TLCondition()
    : ConditionEvaluation(), m_gtAXOL1TLTemplate{nullptr}, m_gtGTB{nullptr}, m_model_wrapper{} {}

l1t::AXOL1TLCondition::AXOL1TLCondition(const GlobalCondition* axol1tlTemplate, const GlobalBoard* ptrGTB)
    : ConditionEvaluation(),
      m_gtAXOL1TLTemplate(static_cast<const AXOL1TLTemplate*>(axol1tlTemplate)),
      m_gtGTB(ptrGTB),
      m_model_wrapper{kModelNamePrefix + m_gtAXOL1TLTemplate->modelVersion()} {}

// copy constructor
void l1t::AXOL1TLCondition::copy(const l1t::AXOL1TLCondition& cp) {
  m_gtAXOL1TLTemplate = cp.gtAXOL1TLTemplate();
  m_gtGTB = cp.gtGTB();

  m_condMaxNumberObjects = cp.condMaxNumberObjects();
  m_condLastResult = cp.condLastResult();
  m_combinationsInCond = cp.getCombinationsInCond();

  m_verbosity = cp.m_verbosity;
}

l1t::AXOL1TLCondition::AXOL1TLCondition(const l1t::AXOL1TLCondition& cp) : ConditionEvaluation() { copy(cp); }

// destructor
l1t::AXOL1TLCondition::~AXOL1TLCondition() {
  // empty
}

// equal operator
l1t::AXOL1TLCondition& l1t::AXOL1TLCondition::operator=(const l1t::AXOL1TLCondition& cp) {
  copy(cp);
  return *this;
}

// methods
void l1t::AXOL1TLCondition::setGtAXOL1TLTemplate(const AXOL1TLTemplate* caloTempl) { m_gtAXOL1TLTemplate = caloTempl; }

///   set the pointer to uGT GlobalBoard
void l1t::AXOL1TLCondition::setuGtB(const GlobalBoard* ptrGTB) { m_gtGTB = ptrGTB; }

/// set score for score saving
void l1t::AXOL1TLCondition::setScore(const float scoreval) const { m_savedscore = scoreval; }

const bool l1t::AXOL1TLCondition::evaluateCondition(const int bxEval) const {
  bool condResult = false;
  int useBx = bxEval + m_gtAXOL1TLTemplate->condRelativeBx();

  // //pointers to objects
  const BXVector<const l1t::Muon*>* candMuVec = m_gtGTB->getCandL1Mu();
  const BXVector<const l1t::L1Candidate*>* candJetVec = m_gtGTB->getCandL1Jet();
  const BXVector<const l1t::L1Candidate*>* candEGVec = m_gtGTB->getCandL1EG();
  const BXVector<const l1t::EtSum*>* candEtSumVec = m_gtGTB->getCandL1EtSum();

  const int NMuons = 4;
  const int NJets = 10;
  const int NEgammas = 4;
  //const int NEtSums = 1;

  //number of indices in vector is #objects * 3 for et, eta, phi
  const int MuVecSize = 12;    //NMuons * 3;      //so 12
  const int JVecSize = 30;     //NJets * 3;        //so 30
  const int EGVecSize = 12;    //NEgammas * 3;    //so 12
  const int EtSumVecSize = 3;  //NEtSums * 3;    //so 3

  //total # inputs in vector is (4+10+4+1)*3 = 57
  const int NInputs = 57;

  //types of inputs and outputs
  typedef ap_fixed<18, 13> inputtype;
  typedef ap_ufixed<18, 14> losstype;

  //define zero
  inputtype fillzero = 0.0;

  //AD vector declaration, will fill later
  inputtype ADModelInput[NInputs] = {};

  //initializing vector by type for my sanity
  inputtype MuInput[MuVecSize];
  inputtype JetInput[JVecSize];
  inputtype EgammaInput[EGVecSize];
  inputtype EtSumInput[EtSumVecSize];

  //declare result vectors +score
  // resulttype result;
  losstype loss;
  // pairtype ADModelResult;  //model outputs a pair of the (result vector, loss)
  float score = -1.0;  //not sure what the best default is hm??

  //check number of input objects we actually have (muons, jets etc)
  int NCandMu = candMuVec->size(useBx);
  int NCandJet = candJetVec->size(useBx);
  int NCandEG = candEGVec->size(useBx);
  int NCandEtSum = candEtSumVec->size(useBx);

  //initialize arrays to zero (std::fill(first, last, value);)
  std::fill(EtSumInput, EtSumInput + EtSumVecSize, fillzero);
  std::fill(MuInput, MuInput + MuVecSize, fillzero);
  std::fill(JetInput, JetInput + JVecSize, fillzero);
  std::fill(EgammaInput, EgammaInput + EGVecSize, fillzero);
  std::fill(ADModelInput, ADModelInput + NInputs, fillzero);

  //then fill the object vectors
  //NOTE assume candidates are already sorted by pt
  //loop over EtSums first, easy because there is max 1 of them
  if (NCandEtSum > 0) {  //check if not empty
    for (int iEtSum = 0; iEtSum < NCandEtSum; iEtSum++) {
      if ((candEtSumVec->at(useBx, iEtSum))->getType() == l1t::EtSum::EtSumType::kMissingEt) {
        EtSumInput[0] =
            ((candEtSumVec->at(useBx, iEtSum))->hwPt()) / 2;  //have to do hwPt/2 in order to match original et inputs
        // EtSumInput[1] = (candEtSumVec->at(useBx, iEtSum))->hwEta(); //this one is zero, so leave it zero
        EtSumInput[2] = (candEtSumVec->at(useBx, iEtSum))->hwPhi();
      }
    }
  }

  //next egammas
  if (NCandEG > 0) {  //check if not empty
    for (int iEG = 0; iEG < NCandEG; iEG++) {
      if (iEG < NEgammas) {  //stop if fill the Nobjects we need
        EgammaInput[0 + (3 * iEG)] = ((candEGVec->at(useBx, iEG))->hwPt()) /
                                     2;  //index 0,3,6,9 //have to do hwPt/2 in order to match original et inputs
        EgammaInput[1 + (3 * iEG)] = (candEGVec->at(useBx, iEG))->hwEta();  //index 1,4,7,10
        EgammaInput[2 + (3 * iEG)] = (candEGVec->at(useBx, iEG))->hwPhi();  //index 2,5,8,11
      }
    }
  }

  //next muons
  if (NCandMu > 0) {  //check if not empty
    for (int iMu = 0; iMu < NCandMu; iMu++) {
      if (iMu < NMuons) {  //stop if fill the Nobjects we need
        MuInput[0 + (3 * iMu)] = ((candMuVec->at(useBx, iMu))->hwPt()) /
                                 2;  //index 0,3,6,9 //have to do hwPt/2 in order to match original et inputs
        MuInput[1 + (3 * iMu)] = (candMuVec->at(useBx, iMu))->hwEtaAtVtx();  //index 1,4,7,10
        MuInput[2 + (3 * iMu)] = (candMuVec->at(useBx, iMu))->hwPhiAtVtx();  //index 2,5,8,11
      }
    }
  }

  //next jets
  if (NCandJet > 0) {  //check if not empty
    for (int iJet = 0; iJet < NCandJet; iJet++) {
      if (iJet < NJets) {  //stop if fill the Nobjects we need
        JetInput[0 + (3 * iJet)] = ((candJetVec->at(useBx, iJet))->hwPt()) /
                                   2;  //index 0,3,6,9...27 //have to do hwPt/2 in order to match original et inputs
        JetInput[1 + (3 * iJet)] = (candJetVec->at(useBx, iJet))->hwEta();  //index 1,4,7,10...28
        JetInput[2 + (3 * iJet)] = (candJetVec->at(useBx, iJet))->hwPhi();  //index 2,5,8,11...29
      }
    }
  }

  //now put it all together-> EtSum+EGamma+Muon+Jet into ADModelInput
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

  //now run the inference
  if ((m_model_wrapper.model_name() == "GTADModel_v3") || (m_model_wrapper.model_name() == "GTADModel_v4")) {
    using resulttype = std::array<ap_fixed<10, 7, AP_RND_CONV, AP_SAT>, 8>;
    loss = readResult<inputtype, resulttype, losstype>(m_model_wrapper, ADModelInput);
  } else {
    using resulttype = ap_fixed<18, 14, AP_RND_CONV, AP_SAT>;
    loss = readResult<inputtype, resulttype, losstype>(m_model_wrapper, ADModelInput);
  }

  // result = ADModelResult.first;
  // loss = ADModelResult.second;
  score = ((loss).to_float()) * 16.0;  //scaling to match threshold
  //save score to class variable in case score saving needed
  setScore(score);

  //number of objects/thrsholds to check
  int iCondition = 0;  // number of conditions: there is only one
  int nObjInCond = m_gtAXOL1TLTemplate->nrObjects();

  if (iCondition >= nObjInCond || iCondition < 0) {
    return false;
  }

  const AXOL1TLTemplate::ObjectParameter objPar = (*(m_gtAXOL1TLTemplate->objectParameter()))[iCondition];

  // condGEqVal indicates the operator used for the condition (>=, =): true for >=
  bool condGEqVal = m_gtAXOL1TLTemplate->condGEq();
  bool passCondition = false;

  passCondition = checkCut(objPar.minAXOL1TLThreshold, score, condGEqVal);

  condResult |= passCondition;  //condresult true if passCondition true else it is false

  //return result
  return condResult;
}

void l1t::AXOL1TLCondition::print(std::ostream& myCout) const {
  myCout << "Dummy Print for AXOL1TLCondition" << std::endl;
  m_gtAXOL1TLTemplate->print(myCout);

  ConditionEvaluation::print(myCout);
}
