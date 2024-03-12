// this class header
#include "L1Trigger/L1TGlobal/interface/ADTTemplate.h"

// system include files
#include <iostream>
#include <iomanip>

ADTTemplate::ADTTemplate() : GlobalCondition() { m_condCategory = l1t::CondADT; }

ADTTemplate::ADTTemplate(const std::string& cName) : GlobalCondition(cName) { m_condCategory = l1t::CondADT; }

ADTTemplate::ADTTemplate(const std::string& cName, const l1t::GtConditionType& cType)  //not sure we need cType
    : GlobalCondition(cName, l1t::CondADT, cType) {
  int nObjects = nrObjects();

  if (nObjects > 0) {
    m_objectType.reserve(nObjects);
  }
}

// copy constructor
ADTTemplate::ADTTemplate(const ADTTemplate& cp) : GlobalCondition(cp.m_condName) { copy(cp); }

// destructor
ADTTemplate::~ADTTemplate() {
  // empty now
}

// assign operator
ADTTemplate& ADTTemplate::operator=(const ADTTemplate& cp) {
  copy(cp);
  return *this;
}

// setConditionParameter - set the parameters of the condition
void ADTTemplate::setConditionParameter(const std::vector<ObjectParameter>& objParameter) {
  m_objectParameter = objParameter;
}

void ADTTemplate::print(std::ostream& myCout) const {
  myCout << "\n  ADTTemplate print..." << std::endl;

  GlobalCondition::print(myCout);

  int nObjects = nrObjects();

  for (int i = 0; i < nObjects; i++) {
    myCout << std::endl;
    myCout << "  Template for object " << i << " [ hex ]" << std::endl;
    myCout << "    ADTThreshold   = " << std::hex << m_objectParameter[i].minADTThreshold << std::endl;
  }

  // reset to decimal output
  myCout << std::dec << std::endl;
}

void ADTTemplate::copy(const ADTTemplate& cp) {
  m_condName = cp.condName();
  m_condCategory = cp.condCategory();
  m_condType = cp.condType();
  m_objectType = cp.objectType();  //not needed for ADT
  m_condGEq = cp.condGEq();
  m_condChipNr = cp.condChipNr();
  m_condRelativeBx = cp.condRelativeBx();

  m_objectParameter = *(cp.objectParameter());
}

// output stream operator
std::ostream& operator<<(std::ostream& os, const ADTTemplate& result) {
  result.print(os);
  return os;
}
