#include <algorithm>
#include <memory>
#include <utility>

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1ScoutingObjectBXVecConverterT.h"

using L1ScoutingCaloTowerBXVecConverter = L1ScoutingObjectBXVecConverterT<l1t::CaloTower, l1ScoutingRun3::CaloTower>;
DEFINE_FWK_MODULE(L1ScoutingCaloTowerBXVecConverter);
