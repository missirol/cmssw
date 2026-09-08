#ifndef L1TriggerScouting_Utilities_convertL1TObjectToL1ScoutingObjectT_h
#define L1TriggerScouting_Utilities_convertL1TObjectToL1ScoutingObjectT_h

#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/TypeDemangler.h"
#include "L1TriggerScouting/Utilities/interface/conversion.h"
#include "L1TriggerScouting/Utilities/interface/convertL1TObjectToL1ScoutingObjectT.h"

template <class T1, class T2>
T2 convertL1TObjectToL1ScoutingObjectT(T1 const&) {
  throw cms::Exception("InvalidType") << "no explicit implementation of convertL1TObjectToL1ScoutingObjectT<"
                                      << edm::typeDemangle(typeid(T1).name()) << ", "
                                      << edm::typeDemangle(typeid(T2).name()) << "> exists!";
}

template <>
inline l1ScoutingRun3::CaloTower convertL1TObjectToL1ScoutingObjectT<l1t::CaloTower, l1ScoutingRun3::CaloTower>(
    l1t::CaloTower const& ct) {
  // Useful links.
  //  - L1T Stage-2 CaloTowerUnpacker:
  //    https://github.com/cms-sw/cmssw/blob/09dbb849b53a0d14d9a8f67556472ba9b5edfba1/EventFilter/L1TRawToDigi/plugins/implementations_stage2/CaloTowerUnpacker.cc#L12
  //  - L1S CaloTower unpacker (2026):
  //    https://github.com/cms-sw/cmssw/blob/09dbb849b53a0d14d9a8f67556472ba9b5edfba1/EventFilter/L1ScoutingRawToDigi/plugins/ScCaloTowerRawToDigi.cc#L94
  //
  // Mask values are hard-coded in order to avoid a circular dependency
  // between this subpackage and EventFilter/L1ScoutingRawToDigi
  // (the header files defining bit shifts, bit masks and conversions for
  // the l1ScoutingRun3:: types should arguably be moved
  // from EventFilter/L1ScoutingRawToDigi to this subpackage).
  //
  // No bit masks are applied in the case of hwEta and hwPhi
  // (they are not expected to be needed, and since ct.hwEta() is
  // a int whose value can be negative, applying a simple mask
  // can lead to misinterpreting its value).
  // Instead, an exception is thrown if either hwEta or hwPhi
  // have invalid values once converted to int16_t.
  int16_t const hwEt = ct.hwPt() & 0x1FF;
  int16_t const erBits = ct.hwEtRatio() & 0x7;
  int16_t const miscBits = ct.hwQual() & 0xF;
  int16_t const hwEta = ct.hwEta();
  int16_t const hwPhi = ct.hwPhi();

  if (not l1ScoutingRun3::calol1::validHwEta(hwEta)) {
    throw cms::Exception("InvalidInput") << "invalid value for CaloTower hwEta (" << hwEta << ")!";
  }

  if (not l1ScoutingRun3::calol1::validHwPhi(hwPhi)) {
    throw cms::Exception("InvalidInput") << "invalid value for CaloTower hwPhi (" << hwPhi << ")!";
  }

  return l1ScoutingRun3::CaloTower(hwEt, erBits, miscBits, hwEta, hwPhi);
}

#endif
