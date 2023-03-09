#include <string>
#include <limits>

#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"
#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "DQM/TrackingMonitor/interface/GetLumi.h"
#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"
#include "DataFormats/OnlineMetaData/interface/OnlineLuminosityRecord.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

namespace {
  typedef dqm::reco::DQMStore DQMStore;
  typedef dqm::reco::MonitorElement MonitorElement;

  struct MEbinning {
    int nbins;
    double xmin;
    double xmax;
  };

  struct Histograms {
    MonitorElement* numberOfPixelClustersVsLS = nullptr;
    MonitorElement* numberOfPixelClustersVsLumi = nullptr;
    MonitorElement* lumiVsLS = nullptr;
    MonitorElement* puVsLS = nullptr;
    MonitorElement* pixelLumiVsLS = nullptr;
    MonitorElement* pixelLumiVsLumi = nullptr;
  };
}  // namespace

//
// class declaration
//

class LumiMonitor : public DQMGlobalEDAnalyzer<Histograms> {
public:
  LumiMonitor(const edm::ParameterSet&);
  ~LumiMonitor() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static void fillHistoPSetDescription(edm::ParameterSetDescription& pset,
                                       int const nbins,
                                       double const xmin,
                                       double const xmax);
  static void fillHistoLSPSetDescription(edm::ParameterSetDescription& pset, int const nbins);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&, Histograms&) const override;
  void dqmAnalyze(edm::Event const&, edm::EventSetup const&, Histograms const&) const override;

  static MEbinning getHistoPSet(const edm::ParameterSet& pset);
  static MEbinning getHistoLSPSet(const edm::ParameterSet& pset);

  std::string const folderName_;

  edm::EDGetTokenT<LumiScalersCollection> const lumiScalersToken_;
  edm::EDGetTokenT<OnlineLuminosityRecord> const onlineMetaDataDigisToken_;
  MEbinning const lumi_binning_;
  MEbinning const ls_binning_;

  bool const doPixelLumi_;
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> const pixelClustersToken_;
  edm::ESGetToken<LHCInfo, LHCInfoRcd> const lhcInfoToken_;
  bool const useBPixLayer1_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> const trkTopoToken_;
  int const minNumberOfPixelsPerCluster_;
  float const minPixelClusterCharge_;
  float const lumi_factor_per_bx_;
};

// -----------------------------
//  constructors and destructor
// -----------------------------

LumiMonitor::LumiMonitor(const edm::ParameterSet& config)
    : folderName_(config.getParameter<std::string>("folderName")),
      lumiScalersToken_(consumes(config.getParameter<edm::InputTag>("scalers"))),
      onlineMetaDataDigisToken_(consumes(config.getParameter<edm::InputTag>("onlineMetaDataDigis"))),
      lumi_binning_(getHistoPSet(
          config.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("lumiPSet"))),
      ls_binning_(getHistoLSPSet(
          config.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("lsPSet"))),
      doPixelLumi_(config.getParameter<bool>("doPixelLumi")),
      pixelClustersToken_(doPixelLumi_ ? consumes<edmNew::DetSetVector<SiPixelCluster>>(
                                             config.getParameter<edm::InputTag>("pixelClusters"))
                                       : edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>>()),
      lhcInfoToken_(doPixelLumi_ ? esConsumes<LHCInfo, LHCInfoRcd>() : edm::ESGetToken<LHCInfo, LHCInfoRcd>()),
      useBPixLayer1_(doPixelLumi_ ? config.getParameter<bool>("useBPixLayer1") : false),
      trkTopoToken_(doPixelLumi_ and (not useBPixLayer1_) ? esConsumes<TrackerTopology, TrackerTopologyRcd>()
                                                          : edm::ESGetToken<TrackerTopology, TrackerTopologyRcd>()),
      minNumberOfPixelsPerCluster_(doPixelLumi_ ? config.getParameter<int>("minNumberOfPixelsPerCluster") : -1),
      minPixelClusterCharge_(doPixelLumi_ ? config.getParameter<double>("minPixelClusterCharge") : -1.),
      lumi_factor_per_bx_(useBPixLayer1_
                              ? GetLumi::FREQ_ORBIT / GetLumi::CM2_TO_NANOBARN / GetLumi::XSEC_PIXEL_CLUSTER
                              : GetLumi::FREQ_ORBIT / GetLumi::CM2_TO_NANOBARN / GetLumi::rXSEC_PIXEL_CLUSTER) {}

MEbinning LumiMonitor::getHistoPSet(const edm::ParameterSet& pset) {
  return MEbinning{
      pset.getParameter<int32_t>("nbins"),
      pset.getParameter<double>("xmin"),
      pset.getParameter<double>("xmax"),
  };
}

MEbinning LumiMonitor::getHistoLSPSet(const edm::ParameterSet& pset) {
  return MEbinning{pset.getParameter<int32_t>("nbins"), -0.5, pset.getParameter<int32_t>("nbins") - 0.5};
}

void LumiMonitor::bookHistograms(DQMStore::IBooker& booker,
                                 edm::Run const& run,
                                 edm::EventSetup const& setup,
                                 Histograms& histograms) const {
  booker.setCurrentFolder(folderName_);

  float const prof_ymin = 0.f;
  float const prof_ymax = std::numeric_limits<float>::max();

  if (doPixelLumi_) {
    auto me = booker.bookProfile("numberOfPixelClustersVsLS",
                                 "number of pixel clusters vs lumisection",
                                 ls_binning_.nbins,
                                 ls_binning_.xmin,
                                 ls_binning_.xmax,
                                 prof_ymin,
                                 prof_ymax);
    me->setAxisTitle("lumisection", 1);
    me->setAxisTitle("number of pixel clusters", 2);
    histograms.numberOfPixelClustersVsLS = me;

    me = booker.bookProfile("numberOfPixelClustersVsLumi",
                            "number of pixel clusters vs online lumi",
                            lumi_binning_.nbins,
                            lumi_binning_.xmin,
                            lumi_binning_.xmax,
                            prof_ymin,
                            prof_ymax);
    me->setAxisTitle("online inst lumi E30 [Hz cm^{-2}]", 1);
    me->setAxisTitle("number of pixel clusters", 2);
    histograms.numberOfPixelClustersVsLumi = me;

    me = booker.bookProfile("pixelLumiVsLS",
                            "pixel-lumi vs lumisection",
                            ls_binning_.nbins,
                            ls_binning_.xmin,
                            ls_binning_.xmax,
                            prof_ymin,
                            prof_ymax);
    me->setAxisTitle("lumisection", 1);
    me->setAxisTitle("pixel-based inst lumi E30 [Hz cm^{-2}]", 2);
    histograms.pixelLumiVsLS = me;

    me = booker.bookProfile("pixelLumiVsLumi",
                            "pixel-lumi vs online lumi",
                            lumi_binning_.nbins,
                            lumi_binning_.xmin,
                            lumi_binning_.xmax,
                            prof_ymin,
                            prof_ymax);
    me->setAxisTitle("online inst lumi E30 [Hz cm^{-2}]", 1);
    me->setAxisTitle("pixel-based inst lumi E30 [Hz cm^{-2}]", 2);
    histograms.pixelLumiVsLumi = me;
  }

  auto me = booker.bookProfile("lumiVsLS",
                               "online lumi vs lumisection",
                               ls_binning_.nbins,
                               ls_binning_.xmin,
                               ls_binning_.xmax,
                               prof_ymin,
                               prof_ymax);
  me->setAxisTitle("lumisection", 1);
  me->setAxisTitle("online inst lumi E30 [Hz cm^{-2}]", 2);
  histograms.lumiVsLS = me;

  me = booker.bookProfile("puVsLS",
                          "online pileup vs lumisection",
                          ls_binning_.nbins,
                          ls_binning_.xmin,
                          ls_binning_.xmax,
                          prof_ymin,
                          prof_ymax);
  me->setAxisTitle("lumisection", 1);
  me->setAxisTitle("online pileup", 2);
  histograms.puVsLS = me;
}

void LumiMonitor::dqmAnalyze(edm::Event const& event,
                             edm::EventSetup const& setup,
                             Histograms const& histograms) const {
  int const ls = event.id().luminosityBlock();

  float online_lumi = -1.f;
  float online_pu = -1.f;
  auto const lumiScalersHandle = event.getHandle(lumiScalersToken_);
  auto const onlineMetaDataDigisHandle = event.getHandle(onlineMetaDataDigisToken_);
  if (lumiScalersHandle.isValid() and not lumiScalersHandle->empty()) {
    auto const scalit = lumiScalersHandle->begin();
    online_lumi = scalit->instantLumi();
    online_pu = scalit->pileup();
  } else if (onlineMetaDataDigisHandle.isValid()) {
    online_lumi = onlineMetaDataDigisHandle->instLumi();
    online_pu = onlineMetaDataDigisHandle->avgPileUp();
  }
  histograms.lumiVsLS->Fill(ls, online_lumi);
  histograms.puVsLS->Fill(ls, online_pu);

  if (doPixelLumi_) {
    size_t pixel_clusters = 0;
    float pixel_lumi = -1.f;
    auto const pixelClusters = event.getHandle(pixelClustersToken_);
    if (pixelClusters.isValid()) {
      TrackerTopology const* tTopo = useBPixLayer1_ ? nullptr : &(setup.getData(trkTopoToken_));
      // Count the number of clusters with at least a minimum
      // number of pixels per cluster and at least a minimum charge.
      size_t tot = 0;
      for (auto pixCluDet = pixelClusters->begin(); pixCluDet != pixelClusters->end(); ++pixCluDet) {
        if (not useBPixLayer1_) {
          DetId const detid = pixCluDet->detId();
          size_t const subdetid = detid.subdetId();
          if (subdetid == PixelSubdetector::PixelBarrel and tTopo->layer(detid) == 1) {
            continue;
          }
        }

        for (auto pixClu = pixCluDet->begin(); pixClu != pixCluDet->end(); ++pixClu) {
          ++tot;
          if (pixClu->size() >= minNumberOfPixelsPerCluster_ and pixClu->charge() >= minPixelClusterCharge_) {
            ++pixel_clusters;
          }
        }
      }

      auto const& lhcInfo = setup.getData(lhcInfoToken_);
edm::LogPrint("") << "AA " << lhcInfo.collidingBunches();
      pixel_lumi = lumi_factor_per_bx_ * lhcInfo.collidingBunches() * pixel_clusters;
    } else {
      pixel_lumi = -1.;
    }

    histograms.numberOfPixelClustersVsLS->Fill(ls, pixel_clusters);
    histograms.numberOfPixelClustersVsLumi->Fill(online_lumi, pixel_clusters);
    histograms.pixelLumiVsLS->Fill(ls, pixel_lumi);
    histograms.pixelLumiVsLumi->Fill(online_lumi, pixel_lumi);
  }
}

void LumiMonitor::fillHistoPSetDescription(edm::ParameterSetDescription& pset,
                                           int const nbins,
                                           double const xmin,
                                           double const xmax) {
  pset.add<int>("nbins", nbins);
  pset.add<double>("xmin", xmin);
  pset.add<double>("xmax", xmax);
}

void LumiMonitor::fillHistoLSPSetDescription(edm::ParameterSetDescription& pset, int const nbins) {
  pset.add<int>("nbins", nbins);
}

void LumiMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pixelClusters", edm::InputTag("hltSiPixelClusters"));
  desc.add<edm::InputTag>("scalers", edm::InputTag("hltScalersRawToDigi"));
  desc.add<edm::InputTag>("onlineMetaDataDigis", edm::InputTag("hltOnlineMetaDataDigis"));
  desc.add<std::string>("folderName", "HLT/LumiMonitoring");
  desc.add<bool>("doPixelLumi", false);
  desc.add<bool>("useBPixLayer1", false);
  desc.add<int>("minNumberOfPixelsPerCluster", 2);  // from DQM/PixelLumi/python/PixelLumiDQM_cfi.py
  desc.add<double>("minPixelClusterCharge", 15000.);

  edm::ParameterSetDescription histoPSet;

  edm::ParameterSetDescription lumiPSet;
  fillHistoPSetDescription(lumiPSet, 9e3, 0, 9e4);  // min/max units: E30 Hz cm^{-2}
  histoPSet.add<edm::ParameterSetDescription>("lumiPSet", lumiPSet);

  edm::ParameterSetDescription lsPSet;
  fillHistoLSPSetDescription(lsPSet, 2500);
  histoPSet.add<edm::ParameterSetDescription>("lsPSet", lsPSet);

  desc.add<edm::ParameterSetDescription>("histoPSet", histoPSet);

  descriptions.addWithDefaultLabel(desc);
}

// Define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LumiMonitor);
