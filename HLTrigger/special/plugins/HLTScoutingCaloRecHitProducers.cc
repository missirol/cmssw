#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"

#include "HLTrigger/special/plugins/HLTScoutingCaloRecHitProducerT.h"

using HLTScoutingEcalRecHitProducer = HLTScoutingCaloRecHitProducerT<EcalRecHitCollection>;
DEFINE_FWK_MODULE(HLTScoutingEcalRecHitProducer);

using HLTScoutingHBHERecHitProducer = HLTScoutingCaloRecHitProducerT<HBHERecHitCollection>;
DEFINE_FWK_MODULE(HLTScoutingHBHERecHitProducer);

using HLTScoutingHORecHitProducer = HLTScoutingCaloRecHitProducerT<HORecHitCollection>;
DEFINE_FWK_MODULE(HLTScoutingHORecHitProducer);

using HLTScoutingHFRecHitProducer = HLTScoutingCaloRecHitProducerT<HFRecHitCollection>;
DEFINE_FWK_MODULE(HLTScoutingHFRecHitProducer);
