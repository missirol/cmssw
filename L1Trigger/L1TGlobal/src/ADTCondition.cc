/**
 * \class ADTCondition
 *
 *
 * Description: evaluation of a condition for adt anomaly detection algorithm
 *
 * Author: Melissa Quinnan
 *
 **/

// this class header
#include "L1Trigger/L1TGlobal/interface/CorrCondition.h"

// system include files
#include <iostream>
#include <iomanip>

#include <string>
#include <vector>
#include <algorithm>
#include "ap_fixed.h"
#include "hls4ml/emulator.h"

// user include files
//   base classes
#include "L1Trigger/L1TGlobal/interface/ADTTemplate.h"
#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"

#include "L1Trigger/L1TGlobal/interface/MuCondition.h"
#include "L1Trigger/L1TGlobal/interface/ADTCondition.h"
#include "L1Trigger/L1TGlobal/interface/CaloCondition.h"
#include "L1Trigger/L1TGlobal/interface/EnergySumCondition.h"
#include "L1Trigger/L1TGlobal/interface/MuonTemplate.h"
#include "L1Trigger/L1TGlobal/interface/CaloTemplate.h"
#include "L1Trigger/L1TGlobal/interface/EnergySumTemplate.h"
#include "L1Trigger/L1TGlobal/interface/GlobalScales.h"

#include "DataFormats/L1Trigger/interface/L1Candidate.h"

#include "L1Trigger/L1TGlobal/interface/GlobalBoard.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

// constructors
//     default
l1t::ADTCondition::ADTCondition() : ConditionEvaluation() {
  // empty
}

//     from base template condition (from event setup usually)
l1t::ADTCondition::ADTCondition(const GlobalCondition* adtTemplate, const GlobalBoard* ptrGTB)
    : ConditionEvaluation(), m_gtADTTemplate(static_cast<const ADTTemplate*>(adtTemplate)), m_gtGTB(ptrGTB) {}

// copy constructor
void l1t::ADTCondition::copy(const l1t::ADTCondition& cp) {
  m_gtADTTemplate = cp.gtADTTemplate();
  m_gtGTB = cp.gtGTB();

  m_condMaxNumberObjects = cp.condMaxNumberObjects();
  m_condLastResult = cp.condLastResult();
  m_combinationsInCond = cp.getCombinationsInCond();

  m_verbosity = cp.m_verbosity;
}

l1t::ADTCondition::ADTCondition(const l1t::ADTCondition& cp) : ConditionEvaluation() { copy(cp); }

// destructor
l1t::ADTCondition::~ADTCondition() {
  // empty
}

// equal operator
l1t::ADTCondition& l1t::ADTCondition::operator=(const l1t::ADTCondition& cp) {
  copy(cp);
  return *this;
}

// methods
void l1t::ADTCondition::setGtADTTemplate(const ADTTemplate* caloTempl) { m_gtADTTemplate = caloTempl; }

///   set the pointer to uGT GlobalBoard
void l1t::ADTCondition::setuGtB(const GlobalBoard* ptrGTB) { m_gtGTB = ptrGTB; }

const bool l1t::ADTCondition::evaluateCondition(const int bxEval) const {
  bool condResult = false;
  int useBx = bxEval + m_gtADTTemplate->condRelativeBx();

  hls4mlEmulator::ModelLoader loader(m_ADTmodelversion);
  std::shared_ptr<hls4mlEmulator::Model> model;

  try {
    model = loader.load_model();
  } catch (std::runtime_error& e) {
    throw cms::Exception("ADTCondition") << "ERROR: failed to load model version " << m_ADTmodelversion;
  }

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
  typedef std::array<ap_fixed<10, 7, AP_RND_CONV, AP_SAT>, 8> resulttype;  //v3
  typedef ap_ufixed<18, 14> losstype;
  typedef std::pair<resulttype, losstype> pairtype;
  // typedef std::array<ap_fixed<10, 7>, 13> resulttype;  //deprecated v1 type:

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
  resulttype result;
  losstype loss;
  pairtype ADModelResult;  //model outputs a pair of the (result vector, loss)
  float score = -1.0;      //not sure what the best default is hm??

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
        MuInput[1 + (3 * iMu)] = (candMuVec->at(useBx, iMu))->hwEta();  //index 1,4,7,10
        MuInput[2 + (3 * iMu)] = (candMuVec->at(useBx, iMu))->hwPhi();  //index 2,5,8,11
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
  model->prepare_input(ADModelInput);  //scaling internal here
  model->predict();
  model->read_result(&ADModelResult);  // this should be the square sum model result

  result = ADModelResult.first;
  loss = ADModelResult.second;
  score = ((loss).to_float()) * 16.0;  //scaling to match threshold

  //number of objects/thrsholds to check
  int iCondition = 0;  // number of conditions: there is only one
  int nObjInCond = m_gtADTTemplate->nrObjects();

  if (iCondition >= nObjInCond || iCondition < 0) {
    return false;
  }

  const ADTTemplate::ObjectParameter objPar = (*(m_gtADTTemplate->objectParameter()))[iCondition];

  // condGEqVal indicates the operator used for the condition (>=, =): true for >=
  bool condGEqVal = m_gtADTTemplate->condGEq();
  bool passCondition = false;

  passCondition = checkCut(objPar.minADTThreshold, score, condGEqVal);

  condResult |= passCondition;  //condresult true if passCondition true else it is false

  //return result
  return condResult;
}

//in order to set model version from config
void l1t::ADTCondition::setModelVersion(const std::string modelversionname) { m_ADTmodelversion = modelversionname; }

void l1t::ADTCondition::print(std::ostream& myCout) const {
  myCout << "Dummy Print for ADTCondition" << std::endl;
  m_gtADTTemplate->print(myCout);

  ConditionEvaluation::print(myCout);
}
