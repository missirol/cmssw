#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "RecoParticleFlow/PFRecHitProducer/interface/JobConfigurationAlpakaRecord.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/PFRecHitHBHETopologyAlpakaESRcd.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/alpaka/PFRecHitHBHETopologyAlpakaESData.h"

//#include "RecoParticleFlow/PFClusterProducer/interface/PFHCALDenseIdNavigatorCore.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFHCALDenseIdNavigator.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFRecHitNavigatorBase.h"

#include <memory>
#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PFRecHitHBHETopologyESProducer : public ESProducer {
  public:
    PFRecHitHBHETopologyESProducer(edm::ParameterSet const& iConfig)
      : hcalEnums_(iConfig.getParameter<std::vector<int>>("hcalEnums")) {
      auto cc = setWhatProduced(this);
      hcalToken_ = cc.consumes();
      geomToken_ = cc.consumes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<std::string>("appendToDataLabel", "");
      desc.add<std::vector<int>>("hcalEnums", {1, 2});
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<PFRecHitHBHETopologyAlpakaESDataHost> produce(PFRecHitHBHETopologyAlpakaESRcd const& iRecord) {
      uint32_t const productSize = 1; //!!

      auto product = std::make_unique<PFRecHitHBHETopologyAlpakaESDataHost>(productSize, cms::alpakatools::host());

      auto const& geom = iRecord.get(geomToken_);
      auto const& topo = iRecord.get(hcalToken_);

      fillPFRecHitHBHETopologyAlpakaESDataHost(product->view(), geom, topo);

      return product;
    }

  private:
    void fillPFRecHitHBHETopologyAlpakaESDataHost(
      PFRecHitHBHETopologyAlpakaESDataHost::View view, CaloGeometry const& geom, HcalTopology const& topo) const;

    std::vector<int> const hcalEnums_;
    edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> hcalToken_;
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  };

  void PFRecHitHBHETopologyESProducer::fillPFRecHitHBHETopologyAlpakaESDataHost(
    PFRecHitHBHETopologyAlpakaESDataHost::View view, CaloGeometry const& geom, HcalTopology const& topo) const {

//    const CaloSubdetectorGeometry* hcalBarrelGeo = geom.getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
//    const CaloSubdetectorGeometry* hcalEndcapGeo = geom.getSubdetectorGeometry(DetId::Hcal, HcalEndcap);
//  
//    // Utilize PFHCALDenseIdNavigatorCore
//    std::unique_ptr<PFHCALDenseIdNavigatorCore> navicore =
//        std::make_unique<PFHCALDenseIdNavigatorCore>(vhcalEnum_, geom, topo);
//  
//    // Filling HCAL DenseID vectors
//    const std::vector<uint32_t> denseId = navicore.get()->getValidDenseIds();

//    // Filling information to define arrays for all relevant HBHE DetIds
//    denseIdMax_ = *max_element(denseId.begin(), denseId.end());
//    denseIdMin_ = *min_element(denseId.begin(), denseId.end());
//    const int detIdArraySize = denseIdMax_ - denseIdMin_ + 1;
//  
//    // Filling detId, positions, neighbours in arrays indexed based on denseId
//    std::vector<uint32_t> detId;
//    detId.clear();
//    detId.resize(detIdArraySize);
//    std::vector<float3> position;
//    position.clear();
//    position.resize(detIdArraySize);
//    std::vector<int> neighbours;
//    neighbours.clear();
//    neighbours.resize(detIdArraySize * 8);
//  
//    for (auto denseid : denseId) {
//      DetId detid = topo.denseId2detId(denseid);
//      HcalDetId hid = HcalDetId(detid);
//      GlobalPoint pos;
//      if (hid.subdet() == HcalBarrel)
//        pos = hcalBarrelGeo->getGeometry(detid)->getPosition();
//      else if (hid.subdet() == HcalEndcap)
//        pos = hcalEndcapGeo->getGeometry(detid)->getPosition();
//      else
//        std::cout << "Unexpected subdetector found for detId " << hid.rawId() << ": " << hid.subdet() << std::endl;
//  
//      unsigned index = getIdx(denseid);
//      detId[index] = (uint32_t)detid;
//      position[index] = make_float3(pos.x(), pos.y(), pos.z());
//  
//      auto neigh = navicore.get()->getNeighbours(denseid);
//  
//      for (uint32_t n = 0; n < 8; n++) {
//        // cmssdt.cern.ch/lxr/source/RecoParticleFlow/PFClusterProducer/interface/PFHCALDenseIdNavigator.h#0087
//        // Order: CENTER(NONE),SOUTH,SOUTHEAST,SOUTHWEST,EAST,WEST,NORTHEAST,NORTHWEST,NORTH
//        // neigh[0] is the rechit itself. Skip for neighbour array
//        // If no neighbour exists in a direction, the value will be 0
//        // Some neighbors from HF included! Need to test if these are included in the map!
//        auto neighDetId = neigh[n + 1].rawId();
//        if (neighDetId > 0 && (&topo)->detId2denseId(neighDetId) >= denseIdMin_ &&
//            (&topo)->detId2denseId(neighDetId) <= denseIdMax_) {
//          neighbours[index * 8 + n] = getIdx(topo.detId2denseId(neighDetId));
//        } else
//          neighbours[index * 8 + n] = -1;
//      }
//    }
//  
//    //
//    // Fill variables for HostAllocator
//    denseId_.resize(denseId.size());
//    std::copy(denseId.begin(), denseId.end(), denseId_.begin());
//    //
//    detId_.resize(detId.size());
//    std::copy(detId.begin(), detId.end(), detId_.begin());
//    neighbours_.resize(neighbours.size());
//    std::copy(neighbours.begin(), neighbours.end(), neighbours_.begin());
//    position_.resize(position.size());
//    std::copy(position.begin(), position.end(), position_.begin());
//
//    navicore.release();
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(PFRecHitHBHETopologyESProducer);
