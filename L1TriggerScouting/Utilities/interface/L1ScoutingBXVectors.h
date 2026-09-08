#ifndef L1TriggerScouting_NanoAOD_L1ScoutingBXVectors_h
#define L1TriggerScouting_NanoAOD_L1ScoutingBXVectors_h

#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCalo.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloJet.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

namespace l1ScoutingRun3 {

  using MuonBxCollection = BXVector<Muon>;
  using JetBxCollection = BXVector<Jet>;
  using EGammaBxCollection = BXVector<EGamma>;
  using TauBxCollection = BXVector<Tau>;
  using BxSumsBxCollection = BXVector<BxSums>;
  using BMTFStubBxCollection = BXVector<BMTFStub>;
  using CaloTowerBxCollection = BXVector<CaloTower>;
  using CaloJetBxCollection = BXVector<CaloJet>;
}  // namespace l1ScoutingRun3

#endif
