#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

//#define LogTrace(X) std::cout << std::endl //!!

class HLTFiltersDQMonitor : public DQMEDAnalyzer {
public:
  explicit HLTFiltersDQMonitor(const edm::ParameterSet&);
  ~HLTFiltersDQMonitor() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void dqmBeginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override;
  void bookHistograms(DQMStore::IBooker&, edm::Run const& iRun, edm::EventSetup const& iSetup) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  bool skipStreamByName(const std::string& streamName) const;
  bool skipPathMonitorElement(const std::string& pathName) const;
  bool skipModuleByEDMType(const std::string& moduleEDMType) const;
  bool skipModuleByType(const std::string& moduleType) const;

  const std::string folderName_;
  const std::string efficPlotNamePrefix_;
  const edm::InputTag triggerResultsInputTag_;
  bool initFailed_;
  bool skipRun_;

  MonitorElement* meMenu_;
  std::unordered_map<std::string, MonitorElement*> meDatasetMap_;
  std::unordered_map<std::string, MonitorElement*> mePathMap_;

  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

  HLTConfigProvider hltConfigProvider_;
};

HLTFiltersDQMonitor::HLTFiltersDQMonitor(const edm::ParameterSet& iConfig)
    : folderName_(iConfig.getParameter<std::string>("folderName")),
      efficPlotNamePrefix_(iConfig.getParameter<std::string>("efficPlotNamePrefix")),
      triggerResultsInputTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),
      initFailed_(false),
      skipRun_(false),
      meMenu_(nullptr) {
  if (triggerResultsInputTag_.process().empty()) {
    edm::LogError("Input") << "Process name not specified in InputTag argument \"triggerResults\""
                           << " (plugin will not produce outputs): \"" << triggerResultsInputTag_.encode() << "\"";
    initFailed_ = true;
    return;
  } else {
    triggerResultsToken_ = consumes<edm::TriggerResults>(triggerResultsInputTag_);
  }
}

void HLTFiltersDQMonitor::dqmBeginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  if (initFailed_) {
    return;
  }

  LogTrace("")
      << "[HLTFiltersDQMonitor] "
      << "----------------------------------------------------------------------------------------------------";
  LogTrace("") << "[HLTFiltersDQMonitor::dqmBeginRun] Run = " << iRun.id();

  // reset data members holding information from the previous run
  skipRun_ = false;

  bool hltChanged(true);
  if (hltConfigProvider_.init(iRun, iSetup, triggerResultsInputTag_.process(), hltChanged)) {
    LogTrace("") << "[HLTFiltersDQMonitor::dqmBeginRun] HLTConfigProvider initialized [processName() = \""
                 << hltConfigProvider_.processName() << "\", tableName() = \"" << hltConfigProvider_.tableName()
                 << "\", size() = " << hltConfigProvider_.size() << "]";
  } else {
    edm::LogError("Input") << "Initialization of HLTConfigProvider failed for Run=" << iRun.id() << " (process=\""
                           << triggerResultsInputTag_.process()
                           << "\") -> plugin will not produce DQM outputs for this run";
    skipRun_ = true;
    return;
  }
}

void HLTFiltersDQMonitor::bookHistograms(DQMStore::IBooker& iBooker,
                                         edm::Run const& iRun,
                                         edm::EventSetup const& iSetup) {
  if (skipRun_ or initFailed_) {
    return;
  }

  iBooker.setCurrentFolder(folderName_);

  iBooker.bookString("HLTMenu", hltConfigProvider_.tableName().c_str());

  auto hltMenuName(hltConfigProvider_.tableName());
  std::replace(hltMenuName.begin(), hltMenuName.end(), '/', '_');
  std::replace(hltMenuName.begin(), hltMenuName.end(), '.', '_');
  while (hltMenuName.front() == '_') {
    hltMenuName.erase(0, 1);
  }

  auto const& triggerNames(hltConfigProvider_.triggerNames());

  meMenu_ = iBooker.bookProfile(efficPlotNamePrefix_ + hltMenuName,
                                "Path Efficiency",
                                triggerNames.size(),
                                0.,
                                triggerNames.size(),
                                -0.1,
                                1.1,
                                "");
  if (meMenu_ and meMenu_->getTProfile() and meMenu_->getTProfile()->GetXaxis()) {
    for (size_t idx = 0; idx < triggerNames.size(); ++idx) {
      meMenu_->getTProfile()->GetXaxis()->SetBinLabel(idx + 1, triggerNames.at(idx).c_str());
    }
  }

  LogTrace("") << "[HLTFiltersDQMonitor::bookHistograms] HLTConfigProvider::size() = " << hltConfigProvider_.size()
               << ", HLTConfigProvider::triggerNames().size() = " << triggerNames.size();

  for (auto const& istream : hltConfigProvider_.streamNames()) {
    LogTrace("") << "[HLTFiltersDQMonitor::bookHistograms] Stream = \"" << istream << "\"";
    if (this->skipStreamByName(istream)) {
      continue;
    }

    auto const& dsets(hltConfigProvider_.streamContent(istream));
    for (auto const& idset : dsets) {
      LogTrace("") << "[HLTFiltersDQMonitor::bookHistograms]   Dataset = \"" << idset << "\"";
      iBooker.setCurrentFolder(folderName_ + "/" + idset);

      const std::vector<std::string>& dsetPathNames = hltConfigProvider_.datasetContent(idset);
      const std::string meDatasetName(efficPlotNamePrefix_ + idset);
      meDatasetMap_[meDatasetName] = iBooker.bookProfile(
          meDatasetName.c_str(), meDatasetName.c_str(), dsetPathNames.size(), 0., dsetPathNames.size(), -0.1, 1.1, "");
      TProfile* meDatasetTProf(nullptr);
      if (meDatasetMap_.at(meDatasetName)) {
        meDatasetTProf = meDatasetMap_.at(meDatasetName)->getTProfile();
      }
      for (size_t idxPath = 0; idxPath < dsetPathNames.size(); ++idxPath) {
        auto const& iPathName(dsetPathNames.at(idxPath));
        if (std::find(triggerNames.begin(), triggerNames.end(), iPathName) == triggerNames.end()) {
          continue;
        }
        if (meDatasetTProf and meDatasetTProf->GetXaxis()) {
          meDatasetTProf->GetXaxis()->SetBinLabel(idxPath + 1, iPathName.c_str());
        }
        if (this->skipPathMonitorElement(iPathName)) {
          continue;
        }

        LogTrace("") << "[HLTFiltersDQMonitor::bookHistograms]     Path = \"" << iPathName << "\"";

        auto const& moduleLabels(hltConfigProvider_.moduleLabels(iPathName));
        std::vector<std::string> mePath_binLabels;
        mePath_binLabels.reserve(moduleLabels.size());
        for (size_t iMod = 0; iMod < moduleLabels.size(); ++iMod) {
          auto const& moduleLabel(moduleLabels.at(iMod));
          if (this->skipModuleByEDMType(hltConfigProvider_.moduleEDMType(moduleLabel)) or
              this->skipModuleByType(hltConfigProvider_.moduleType(moduleLabel))) {
            LogTrace("") << "[HLTFiltersDQMonitor::bookHistograms]       [-] Module = \"" << moduleLabel << "\"";
            continue;
          }
          LogTrace("") << "[HLTFiltersDQMonitor::bookHistograms]       [bin=" << mePath_binLabels.size() + 1
                       << "] Module = \"" << moduleLabel << "\"";
          mePath_binLabels.emplace_back(moduleLabel);
        }

        if (mePath_binLabels.empty()) {
          continue;
        }

        const std::string mePathName(efficPlotNamePrefix_ + idset + "_" + iPathName);
        mePathMap_[mePathName] = iBooker.bookProfile(
            mePathName.c_str(), iPathName.c_str(), mePath_binLabels.size(), 0., mePath_binLabels.size(), -0.1, 1.1, "");

        if (mePathMap_.at(mePathName)) {
          auto* const mePathTProf(mePathMap_.at(mePathName)->getTProfile());
          if (mePathTProf and mePathTProf->GetXaxis()) {
            for (size_t iMod = 0; iMod < mePath_binLabels.size(); ++iMod) {
              mePathTProf->GetXaxis()->SetBinLabel(iMod + 1, mePath_binLabels.at(iMod).c_str());
            }
          }
        }
      }
    }
  }
}

void HLTFiltersDQMonitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (skipRun_ or initFailed_) {
    return;
  }

  LogTrace("") << "[HLTFiltersDQMonitor::analyze] --------------------------------------------------------";
  LogTrace("") << "[HLTFiltersDQMonitor::analyze] Run = " << iEvent.id().run()
               << ", LuminosityBlock = " << iEvent.id().luminosityBlock() << ", Event = " << iEvent.id().event();

  auto const& triggerResults(iEvent.getHandle(triggerResultsToken_));

  if (not triggerResults.isValid()) {
    edm::LogWarning("Input") << "Invalid handle to edm::TriggerResults (InputTag: \"triggerResults\")"
                             << " -> plugin will not fill DQM outputs for this event";
    return;
  }

  // fill MonitorElement: HLT-Menu (bin: path)
  if (meMenu_ and meMenu_->getTProfile() and meMenu_->getTProfile()->GetXaxis()) {
    auto const& triggerNames(hltConfigProvider_.triggerNames());
    for (auto const& iPathName : triggerNames) {
      const uint pathIndex(hltConfigProvider_.triggerIndex(iPathName));
      if (pathIndex >= triggerResults->size()) {
        edm::LogError("Logic") << "Index associated to path \"" << iPathName << "\" (" << pathIndex
                               << ") is inconsistent with triggerResults::size() (" << triggerResults->size()
                               << ") -> plugin will not fill bin associated to this path in HLT-Menu MonitorElement";
        continue;
      }
      auto const pathAccept(triggerResults->accept(pathIndex));
      auto* const axis(meMenu_->getTProfile()->GetXaxis());
      auto const ibin(axis->FindBin(iPathName.c_str()));
      if ((0 < ibin) and (ibin <= axis->GetNbins())) {
        meMenu_->Fill(ibin - 0.5, pathAccept);
      }
    }
  }

  // fill MonitorElements for PrimaryDatasets and Paths
  // loop over Streams
  for (auto const& istream : hltConfigProvider_.streamNames()) {
    LogTrace("") << "[HLTFiltersDQMonitor::analyze]   Stream = \"" << istream << "\"";
    auto const& dsets(hltConfigProvider_.streamContent(istream));
    // loop over PrimaryDatasets in Stream
    for (auto const& idset : dsets) {
      LogTrace("") << "[HLTFiltersDQMonitor::analyze]     Dataset = \"" << idset << "\"";
      TProfile* meDatasetTProf(nullptr);
      const std::string meDatasetName(efficPlotNamePrefix_ + idset);
      if (meDatasetMap_.find(meDatasetName) != meDatasetMap_.end()) {
        meDatasetTProf = meDatasetMap_.at(meDatasetName)->getTProfile();
      }
      auto const& dsetPathNames(hltConfigProvider_.datasetContent(idset));
      // loop over Paths in PrimaryDataset
      for (auto const& iPathName : dsetPathNames) {
        const uint pathIndex(hltConfigProvider_.triggerIndex(iPathName));
        if (pathIndex >= triggerResults->size()) {
          edm::LogError("Logic") << "Index associated to path \"" << iPathName << "\" (" << pathIndex
                                 << ") is inconsistent with triggerResults::size() (" << triggerResults->size()
                                 << ") -> plugin will not fill DQM info related to this path";
          continue;
        }
        auto const pathAccept(triggerResults->accept(pathIndex));
        LogTrace("") << "[HLTFiltersDQMonitor::analyze]       "
                     << "Path = \"" << iPathName << "\", HLTConfigProvider::triggerIndex(\"" << iPathName
                     << "\") = " << pathIndex << ", Accept = " << pathAccept;
        // fill MonitorElement: PrimaryDataset (bin: path)
        if (meDatasetTProf and meDatasetTProf->GetXaxis()) {
          auto* const axis(meDatasetTProf->GetXaxis());
          auto const ibin(axis->FindBin(iPathName.c_str()));
          if ((0 < ibin) and (ibin <= axis->GetNbins())) {
            meDatasetTProf->Fill(ibin - 0.5, pathAccept);
          }
        }
        // fill MonitorElement: Path (bin: filter)
        auto const mePathName(efficPlotNamePrefix_ + idset + "_" + iPathName);
        if (mePathMap_.find(mePathName) != mePathMap_.end()) {
          auto* const mePathTProf(mePathMap_.at(mePathName)->getTProfile());
          if (mePathTProf) {
            auto* const axis(mePathTProf->GetXaxis());
            if (axis) {
              // index of module giving path decision (last filter executed in the path)
              const auto lastModuleIndex(triggerResults->index(pathIndex));
              const auto lastModuleLabel(hltConfigProvider_.moduleLabel(pathIndex, lastModuleIndex));
              LogTrace("") << "[HLTFiltersDQMonitor::analyze]         "
                           << "lastModuleLabel = " << lastModuleLabel;
              // number of modules in the path
              const unsigned sizeModulesPath(hltConfigProvider_.size(pathIndex));
              LogTrace("") << "[HLTFiltersDQMonitor::analyze]         "
                           << "-> selected lastModuleIndex = " << lastModuleIndex << " (HLTConfigProvider::size("
                           << pathIndex << ") = " << sizeModulesPath << ")";
              if (lastModuleIndex >= sizeModulesPath) {
                edm::LogError("Logic") << "Selected index (" << lastModuleIndex << ") for last filter of path \""
                                       << iPathName << "\" is inconsistent with number of modules in the path ("
                                       << sizeModulesPath << ")";
                continue;
              }
              // filter-bins are filled, with a 0 or 1, only when all previous filters in the path have passed
              for (size_t modIdx = 0; modIdx <= lastModuleIndex; ++modIdx) {
                // consider only selected EDFilter modules
                auto const& moduleLabel(hltConfigProvider_.moduleLabel(pathIndex, modIdx));
                if (this->skipModuleByEDMType(hltConfigProvider_.moduleEDMType(moduleLabel)) or
                    this->skipModuleByType(hltConfigProvider_.moduleType(moduleLabel))) {
                  continue;
                }
                auto const filterAccept((modIdx == lastModuleIndex) ? pathAccept : true);
                LogTrace("") << "[HLTFiltersDQMonitor::analyze]         "
                             << "HLTConfigProvider::moduleLabel(" << pathIndex << ", " << modIdx << ") = \""
                             << moduleLabel << "\", filterAccept = " << filterAccept;
                auto const ibin(axis->FindBin(moduleLabel.c_str()));
                if ((0 < ibin) and (ibin <= axis->GetNbins())) {
                  mePathTProf->Fill(ibin - 0.5, filterAccept);
                }
              }
            }
          }
        }
      }
    }
  }
}

bool HLTFiltersDQMonitor::skipStreamByName(const std::string& streamName) const {
  if ((streamName.find("Physics") != std::string::npos) || (streamName.find("Scouting") != std::string::npos) ||
      (streamName.find("Parking") != std::string::npos) || (streamName == "A")) {
    return false;
  }
  return true;
}

bool HLTFiltersDQMonitor::skipPathMonitorElement(const std::string& pathName) const {
  if ((pathName.find("HLT_") == std::string::npos) || (pathName.find("HLT_Physics") != std::string::npos) ||
      (pathName.find("HLT_Random") != std::string::npos)) {
    return true;
  }
  return false;
}

bool HLTFiltersDQMonitor::skipModuleByEDMType(const std::string& moduleEDMType) const {
  return (moduleEDMType != "EDFilter");
}

bool HLTFiltersDQMonitor::skipModuleByType(const std::string& moduleType) const { return (moduleType == "HLTBool"); }

void HLTFiltersDQMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("folderName", "HLT/Filters");
  desc.add<std::string>("efficPlotNamePrefix", "effic_");
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  descriptions.add("hltFiltersDQMonitor", desc);
}

DEFINE_FWK_MODULE(HLTFiltersDQMonitor);
