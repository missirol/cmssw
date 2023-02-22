#ifndef HeterogeneousCore_AlpakaTest_interface_AlpakaESTestData_h
#define HeterogeneousCore_AlpakaTest_interface_AlpakaESTestData_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "RecoParticleFlow/PFRecHitProducer/interface/PFRecHitHBHEParamsAlpakaESDataSoA.h"

namespace reco {

  using PFRecHitHBHEParamsAlpakaESDataHost = PortableHostCollection<PFRecHitHBHEParamsAlpakaESDataSoA>;

}

#endif
