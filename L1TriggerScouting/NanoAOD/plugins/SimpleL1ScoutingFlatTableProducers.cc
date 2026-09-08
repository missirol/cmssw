#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCalo.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloJet.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

using SimpleL1ScoutingBMTFStubFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::BMTFStub>;
DEFINE_FWK_MODULE(SimpleL1ScoutingBMTFStubFlatTableProducer);

using SimpleL1ScoutingBxSumsFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::BxSums>;
DEFINE_FWK_MODULE(SimpleL1ScoutingBxSumsFlatTableProducer);

using SimpleL1ScoutingCaloJetFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::CaloJet>;
DEFINE_FWK_MODULE(SimpleL1ScoutingCaloJetFlatTableProducer);

using SimpleL1ScoutingCaloTowerFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::CaloTower>;
DEFINE_FWK_MODULE(SimpleL1ScoutingCaloTowerFlatTableProducer);

using SimpleL1ScoutingEGammaFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::EGamma>;
DEFINE_FWK_MODULE(SimpleL1ScoutingEGammaFlatTableProducer);

using SimpleL1ScoutingJetFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::Jet>;
DEFINE_FWK_MODULE(SimpleL1ScoutingJetFlatTableProducer);

using SimpleL1ScoutingMuonFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::Muon>;
DEFINE_FWK_MODULE(SimpleL1ScoutingMuonFlatTableProducer);

using SimpleL1ScoutingTauFlatTableProducer = BXVectorSimpleFlatTableProducer<l1ScoutingRun3::Tau>;
DEFINE_FWK_MODULE(SimpleL1ScoutingTauFlatTableProducer);
